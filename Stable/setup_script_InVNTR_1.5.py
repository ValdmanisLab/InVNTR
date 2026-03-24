#!/usr/bin/env python3
import csv
import os
import subprocess
import requests
import shutil
import argparse

# --- Configuration ---
GITHUB_RAW_BASE = "https://raw.githubusercontent.com/ValdmanisLab/InVNTR/main/Stable/setup"
CSV_REGULAR     = "assemblies.csv"
CSV_TESTING     = "assemblies_testing.csv"

# 1000 Genomes config
KG_BASE_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes"
KG_CHRS = [str(i) for i in range(1, 23)] + ["X"]

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
    print(f"⬇️  Downloading {url} ...")
    if url.startswith("s3://"):
        subprocess.run(["aws", "s3", "cp", url, dest], check=False)
    else:
        r = requests.get(url, stream=True)
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

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
    print(f"    🔓 Unzipping {gz_path} with pigz...")
    subprocess.run(["pigz", "-d", "-f", gz_path], check=True)
    return unzipped_name

# --- 1000 Genomes download ---
def download_1000g_vcfs():
    print("\n🧬 Downloading 1000 Genomes GRCh38 phased VCFs...")
    for chr_ in KG_CHRS:
        for suffix in ["vcf.gz", "vcf.gz.tbi"]:
            filename = f"ALL.chr{chr_}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.{suffix}"
            url = f"{KG_BASE_URL}/{filename}"
            dest = filename  # current directory
            if os.path.exists(dest):
                print(f"    ⚠️ {dest} exists, verifying size...")
            download_file(url, dest)
            local_size = file_size(dest)
            remote_size = http_file_size(url)
            if remote_size and local_size == remote_size:
                print(f"    ✅ Verified {filename} (size match)")
            else:
                print(f"    ⚠️ Size mismatch for {filename} (local {local_size}, remote {remote_size})")

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

    if not missing:
        print("✅ All files present.\n")
    else:
        print(f"⚠️ {len(missing)} missing file(s). Downloading...\n")
        for old, url in missing:
            if not url:
                print(f"⚠️ No online location for {old}")
                continue
            download_file(url, old)
            local_size = file_size(old)
            remote_size = s3_file_size(url) if url.startswith("s3://") else http_file_size(url)
            if remote_size and local_size == remote_size:
                print(f"✅ Verified {old} (size match)")
            else:
                print(f"⚠️ Size mismatch for {old} (local {local_size}, remote {remote_size})")

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

    # Step 3: Unzip .gz files with pigz and rename
    print("\n🗜️  Unzipping and renaming .gz files with pigz...")
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