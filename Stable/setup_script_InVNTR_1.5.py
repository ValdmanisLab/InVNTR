#!/usr/bin/env python3
import csv
import os
import subprocess
import requests
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# --- Configuration ---
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/ValdmanisLab/InVNTR/main/Stable/setup"
CSV_REGULAR     = "assemblies.csv"
CSV_TESTING     = "assemblies_testing.csv"

# 1000 Genomes — single sites-only VCF (2 files: VCF + index)
KG_SITES_VCF = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
KG_SITES_TBI = KG_SITES_VCF + ".tbi"

# Max parallel downloads
MAX_WORKERS = 4
MAX_RETRIES = 3

print_lock = threading.Lock()

def tprint(*args, **kwargs):
    """Thread-safe print."""
    with print_lock:
        print(*args, **kwargs)

# --- Argument parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Download and prepare genomic reference data for InVNTR."
    )
    parser.add_argument(
        "-t", "--testing",
        action="store_true",
        help="Use the smaller testing manifest (assemblies_testing.csv) instead of the full one."
    )
    return parser.parse_args()

# --- CSV download ---
def download_csv(testing=False):
    csv_name = CSV_TESTING if testing else CSV_REGULAR
    url = f"{GITHUB_RAW_BASE}/{csv_name}"
    print(f"📥 Downloading manifest: {csv_name} ...")
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(csv_name, "w", encoding="utf-8") as f:
            f.write(r.text)
        print(f"✅ Saved {csv_name} to current directory.\n")
        return csv_name
    except requests.RequestException as e:
        print(f"❌ Failed to download {csv_name}: {e}")
        raise SystemExit(1)

# --- Utility functions ---
def file_size(path):
    try:
        return os.path.getsize(path)
    except FileNotFoundError:
        return 0

def s3_file_size(url):
    try:
        result = subprocess.run(
            ["aws", "s3", "ls", url],
            capture_output=True, text=True, check=True
        )
        parts = result.stdout.split()
        if len(parts) >= 3:
            return int(parts[2])
    except subprocess.CalledProcessError:
        pass
    return None

def http_file_size(url):
    try:
        r = requests.head(url, allow_redirects=True, timeout=10)
        if "Content-Length" in r.headers:
            return int(r.headers["Content-Length"])
    except Exception:
        pass
    return None

def download_file(url, dest):
    tprint(f"⬇️  Downloading {url} ...")
    if url.startswith("s3://"):
        subprocess.run(["aws", "s3", "cp", url, dest], check=False)
    else:
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                r = requests.get(url, stream=True, timeout=60)
                with open(dest, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return  # success
            except (requests.exceptions.ChunkedEncodingError,
                    requests.exceptions.ConnectionError,
                    requests.exceptions.Timeout) as e:
                if attempt < MAX_RETRIES:
                    tprint(f"    ⚠️  Attempt {attempt}/{MAX_RETRIES} failed for {os.path.basename(dest)}: {e}. Retrying...")
                else:
                    tprint(f"    ❌ All {MAX_RETRIES} attempts failed for {os.path.basename(dest)}")
                    raise

def download_and_verify(filename, url):
    """Download a single file and verify its size. Returns (filename, success, message)."""
    if not url:
        return filename, False, f"⚠️  No online location for {filename}"
    try:
        download_file(url, filename)
    except Exception as e:
        return filename, False, f"❌ Download failed for {filename}: {e}"
    local_size = file_size(filename)
    remote_size = s3_file_size(url) if url.startswith("s3://") else http_file_size(url)
    if remote_size and local_size == remote_size:
        return filename, True, f"✅ Verified {filename} (size match)"
    else:
        return filename, False, f"⚠️  Size mismatch for {filename} (local {local_size}, remote {remote_size})"

def generate_fai(fasta_path):
    """Generate .fai index using samtools faidx."""
    if not os.path.exists(fasta_path):
        print(f"⚠️ Skipping {fasta_path} — file not found.")
        return False
    fai_path = fasta_path + ".fai"
    if os.path.exists(fai_path):
        print(f"    ✅ {fai_path} already exists.")
        return True
    try:
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
        print(f"    🧬 Created index: {fai_path}")
        return True
    except subprocess.CalledProcessError:
        print(f"    ⚠️ Failed to create index for {fasta_path}")
        return False

def unzip_with_pigz(gz_path):
    if not os.path.exists(gz_path):
        return None
    unzipped_name = gz_path[:-3]
    # Prefer pigz (faster, parallel), fall back to gzip
    if shutil.which("pigz"):
        tool = ["pigz", "-d", "-f"]
        tool_name = "pigz"
    elif shutil.which("gzip"):
        tool = ["gzip", "-d", "-f"]
        tool_name = "gzip"
    else:
        raise RuntimeError("Neither pigz nor gzip found — please install one (e.g. sudo apt install pigz)")
    print(f"    🔓 Unzipping {gz_path} with {tool_name}...")
    subprocess.run(tool + [gz_path], check=True)
    return unzipped_name

# --- 1000 Genomes download ---
def download_1000g_vcfs():
    print("\n🧬 Downloading 1000 Genomes GRCh38 sites-only VCF...")

    vcf_jobs = []
    for url in [KG_SITES_VCF, KG_SITES_TBI]:
        filename = url.split("/")[-1]
        if os.path.exists(filename):
            print(f"    ✅ {filename} already present, skipping.")
        else:
            vcf_jobs.append((filename, url))

    if not vcf_jobs:
        print("✅ All VCF files already present.")
        return

    print(f"⬇️  Downloading {len(vcf_jobs)} VCF file(s)...\n")

    completed = 0
    with ThreadPoolExecutor(max_workers=min(MAX_WORKERS, len(vcf_jobs))) as executor:
        futures = {executor.submit(download_and_verify, filename, url): filename
                   for filename, url in vcf_jobs}
        for future in as_completed(futures):
            filename, success, msg = future.result()
            completed += 1
            tprint(f"  [{completed}/{len(vcf_jobs)}] {msg}")

# --- Main logic for assemblies ---
def process_assemblies(csv_file):
    print("🔍 Reading CSV and checking files...\n")
    missing = []
    renamed = 0
    unzipped = 0
    indexed = 0

    with open(csv_file, newline='', encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)

    # Step 1: Find missing files
    for row in rows:
        old = row.get("Old_Filename", "").strip()
        new = row.get("Filename", "").strip()
        url = row.get("online_location", "").strip()
        if not old:
            continue
        if not os.path.exists(old) and not os.path.exists(new):
            print(f"❌ Missing {old}")
            missing.append((old, url))

    # Count VCF files that will be downloaded later, for the grand total
    vcf_needed = sum(
        1 for url in [KG_SITES_VCF, KG_SITES_TBI]
        if not os.path.exists(url.split("/")[-1])
    )

    total_to_download = len(missing) + vcf_needed
    print(f"\n📊 Download summary:")
    print(f"   Assembly files missing : {len(missing)}")
    print(f"   VCF files missing      : {vcf_needed}")
    print(f"   Total files to download: {total_to_download}\n")

    if not missing:
        print("✅ All assembly files present.\n")
    else:
        s3_jobs   = [(old, url) for old, url in missing if url.startswith("s3://")]
        http_jobs = [(old, url) for old, url in missing if not url.startswith("s3://") and url]
        no_url    = [(old, url) for old, url in missing if not url]

        for old, _ in no_url:
            print(f"⚠️  No online location for {old}")

        print(f"⬇️  Downloading {len(s3_jobs)} S3 files and {len(http_jobs)} HTTP files "
              f"({MAX_WORKERS} at a time)...\n")

        completed = 0
        total_assembly = len(s3_jobs) + len(http_jobs)
        all_jobs = [(old, url) for old, url in missing if url]

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(download_and_verify, old, url): old
                       for old, url in all_jobs}
            for future in as_completed(futures):
                old, success, msg = future.result()
                completed += 1
                tprint(f"  [{completed}/{total_assembly}] {msg}")

    # Step 2: Rename uncompressed FASTA files
    print("\n🔄 Renaming uncompressed .fasta files...")
    for row in rows:
        old = row.get("Old_Filename", "").strip()
        new = row.get("Filename", "").strip()
        if old.endswith(".gz") or not old:
            continue
        if os.path.exists(old) and not os.path.exists(new):
            os.rename(old, new)
            print(f"    🔁 {old} → {new}")
            renamed += 1
    print(f"✅ Renamed {renamed} file(s).")

    # Step 3: Unzip .gz files and rename
    print("\n🗜️  Unzipping and renaming .gz files...")
    for row in rows:
        old = row.get("Old_Filename", "").strip()
        new = row.get("Filename", "").strip()
        if not old.endswith(".gz"):
            continue
        if os.path.exists(old):
            unzipped_name = unzip_with_pigz(old)
            if unzipped_name:
                os.rename(unzipped_name, new)
                print(f"    🔁 {unzipped_name} → {new}")
                unzipped += 1
    print(f"✅ Unzipped and renamed {unzipped} file(s).")

    # Step 4: Final verification
    print("\n🧾 Verifying final file list...")
    missing_final = 0
    for row in rows:
        new = row.get("Filename", "").strip()
        if not os.path.exists(new):
            print(f"❌ Missing final file: {new}")
            missing_final += 1
    if missing_final == 0:
        print("✅ All files accounted for and properly named.")
    else:
        print(f"⚠️ {missing_final} file(s) still missing at the end.")

    # Step 5: Generate .fai for reference sequences
    print("\n🧬 Generating .fai index files for reference (Haplotype == 'ref')...")
    for row in rows:
        hap = row.get("Haplotype", "").strip().lower()
        fasta = row.get("Filename", "").strip()
        if hap == "ref" and fasta and os.path.exists(fasta):
            if generate_fai(fasta):
                indexed += 1
    print(f"✅ Indexed {indexed} reference file(s).")

# --- Main entry ---
def main():
    args = parse_args()
    csv_file = download_csv(testing=args.testing)
    process_assemblies(csv_file)
    download_1000g_vcfs()

if __name__ == "__main__":
    main()