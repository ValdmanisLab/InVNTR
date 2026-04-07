#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
GITHUB_REPO="ValdmanisLab/InVNTR"
GITHUB_RAW_BASE="https://raw.githubusercontent.com/ValdmanisLab/InVNTR/main/Stable/setup"
CSV_REGULAR="assemblies.csv"
CSV_TESTING="assemblies_testing.csv"

AGC_FULL_URL="https://uwnetid.sharepoint.com/:u:/s/valdmanislabcollab/IQC0v8k0UAcHTIoAp1nyefjzAfQuh_ILvfv5cs-fGuCNA8c?e=Ma0EQo&download=1"
AGC_TESTING_URL="https://uwnetid.sharepoint.com/:u:/s/valdmanislabcollab/IQCpBMmpRvviR6sMnaoUu94OAcQi-DdMtRcXWL9HRH-Nyx8?e=c9gKWl&download=1"
AGC_FULL_NAME="1546.agc"
AGC_TESTING_NAME="1546_testing.agc"

KG_SITES_VCF="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
KG_SITES_TBI="${KG_SITES_VCF}.tbi"

CONDA_ENV_NAME="invntr"
LR_DIR="LR"
TESTING=false

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--testing)
            TESTING=true
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash setup.sh [-t|--testing]"
            exit 1
            ;;
    esac
done

if [ "$TESTING" = true ]; then
    CSV_NAME="$CSV_TESTING"
    AGC_NAME="$AGC_TESTING_NAME"
    AGC_URL="$AGC_TESTING_URL"
else
    CSV_NAME="$CSV_REGULAR"
    AGC_NAME="$AGC_FULL_NAME"
    AGC_URL="$AGC_FULL_URL"
fi

# --- Helpers ---
log()  { echo "  $*"; }
ok()   { echo "  ✅ $*"; }
info() { echo "  ⬇️  $*"; }
warn() { echo "  ⚠️  $*"; }
fail() { echo "  ❌ $*"; exit 1; }

wget_file() {
    local url="$1"
    local dest="$2"
    if [ -f "$dest" ]; then
        ok "$(basename "$dest") already present, skipping."
        return
    fi
    info "Downloading $(basename "$dest") ..."
    wget -O "$dest" "$url" || fail "wget failed for $(basename "$dest")"
    ok "Downloaded $(basename "$dest")."
}

# --- Step 0: Create LR/ directory ---
echo ""
echo "📁 Checking LR/ directory..."
if [ ! -d "$LR_DIR" ]; then
    mkdir -p "$LR_DIR"
    ok "Created $LR_DIR/"
else
    ok "$LR_DIR/ already exists."
fi

# --- Step 1: Find mamba or conda ---
echo ""
echo "🐍 Looking for mamba or conda..."
if command -v mamba &> /dev/null; then
    CONDA_TOOL="mamba"
    ok "Found mamba."
elif command -v conda &> /dev/null; then
    CONDA_TOOL="conda"
    ok "Found conda (mamba not available)."
else
    fail "Neither mamba nor conda found. Please install one before running this script."
fi

# --- Step 2: Download latest release assets (.py and .yml) ---
echo ""
echo "🔽 Fetching latest InVNTR release from GitHub..."
RELEASE_JSON=$(curl -fsSL "https://api.github.com/repos/${GITHUB_REPO}/releases/latest")
RELEASE_TAG=$(echo "$RELEASE_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin)['tag_name'])")
log "Latest release: $RELEASE_TAG"

ENV_YML=""
# Parse asset names and download URLs using python3 for reliable JSON parsing
while IFS=$'\t' read -r name url; do
    if [[ "$name" == *.py ]] || [[ "$name" == *.yml ]]; then
        if [ -f "$name" ]; then
            ok "$name already exists, skipping."
        else
            info "Downloading $name ..."
            wget -O "$name" "$url" || fail "wget failed for $name"
            ok "Downloaded $name."
        fi
        if [[ "$name" == *.yml ]]; then
            ENV_YML="$name"
        fi
    fi
done < <(echo "$RELEASE_JSON" | python3 -c "
import sys, json
assets = json.load(sys.stdin).get('assets', [])
for a in assets:
    print(a['name'] + '\t' + a['browser_download_url'])
")

if [ -z "$ENV_YML" ]; then
    fail "No .yml file found in the latest release. Cannot create conda environment."
fi

# --- Step 3: Create conda environment ---
echo ""
echo "🐍 Checking for '$CONDA_ENV_NAME' conda environment..."
if $CONDA_TOOL env list | grep -q "^${CONDA_ENV_NAME}\s"; then
    ok "'$CONDA_ENV_NAME' environment already exists, skipping creation."
else
    log "Creating '$CONDA_ENV_NAME' environment from $ENV_YML ..."
    $CONDA_TOOL env create -f "$ENV_YML" || fail "Failed to create conda environment."
    ok "Environment '$CONDA_ENV_NAME' created."
fi

# Activate the environment
log "Activating '$CONDA_ENV_NAME' environment..."
# Source conda init so that activation works in a non-interactive shell
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_NAME"
ok "Environment '$CONDA_ENV_NAME' is active."

# --- Step 4: Download CSV manifest ---
echo ""
echo "📥 Downloading manifest: $CSV_NAME ..."
wget_file "${GITHUB_RAW_BASE}/${CSV_NAME}" "${LR_DIR}/${CSV_NAME}"
CSV_PATH="${LR_DIR}/${CSV_NAME}"

# --- Step 5: Check if assemblies already extracted ---
echo ""
echo "📂 Checking if assemblies already present in $LR_DIR/..."

# Read Filename column from CSV (skip header, handle UTF-8 BOM)
FILENAMES=$(python3 -c "
import csv, sys
with open('${CSV_PATH}', newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fn = row.get('Filename','').strip()
        if fn:
            print(fn)
")

ALL_PRESENT=true
while IFS= read -r fname; do
    if [ ! -f "${LR_DIR}/${fname}" ]; then
        ALL_PRESENT=false
        break
    fi
done <<< "$FILENAMES"

if [ "$ALL_PRESENT" = true ]; then
    ok "All assembly files already present in $LR_DIR/, skipping .agc download and extraction."
else
    # --- Step 5a: Download .agc ---
    echo ""
    echo "📦 Downloading .agc archive: $AGC_NAME ..."
    AGC_PATH="${LR_DIR}/${AGC_NAME}"
    wget_file "$AGC_URL" "$AGC_PATH"

    # --- Step 5b: Extract assemblies ---
    echo ""
    echo "🗜️  Extracting assemblies from $AGC_PATH ..."
    SAMPLE_NAMES=$(python3 -c "
import csv, sys
with open('${CSV_PATH}', newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fn = row.get('Filename','').strip()
        if fn:
            name = fn[:-3] if fn.endswith('.fa') else fn
            print(name)
")

    TOTAL=$(echo "$SAMPLE_NAMES" | wc -l | tr -d ' ')
    COUNT=0
    FAILED=()

    while IFS= read -r name; do
        COUNT=$((COUNT + 1))
        OUT_PATH="${LR_DIR}/${name}.fa"
        if [ -f "$OUT_PATH" ]; then
            ok "[$COUNT/$TOTAL] $OUT_PATH already exists, skipping."
            continue
        fi
        log "[$COUNT/$TOTAL] Extracting $name → $OUT_PATH ..."
        if agc getset "$AGC_PATH" "$name" > "$OUT_PATH"; then
            ok "Extracted $name."
        else
            warn "Failed to extract $name."
            FAILED+=("$name")
        fi
    done <<< "$SAMPLE_NAMES"

    if [ ${#FAILED[@]} -gt 0 ]; then
        warn "${#FAILED[@]} sample(s) failed to extract: ${FAILED[*]}"
    else
        ok "Extraction complete — $TOTAL file(s) written."
    fi

    # --- Step 5c: Delete .agc ---
    log "Removing $AGC_PATH ..."
    rm -f "$AGC_PATH" && ok "Deleted $AGC_PATH." || warn "Could not delete $AGC_PATH."
fi

# --- Step 6: Download 1000 Genomes VCF files ---
echo ""
echo "🧬 Downloading 1000 Genomes GRCh38 sites-only VCF..."
VCF_FILENAME=$(basename "$KG_SITES_VCF")
TBI_FILENAME=$(basename "$KG_SITES_TBI")
wget_file "$KG_SITES_VCF" "${LR_DIR}/${VCF_FILENAME}"
wget_file "$KG_SITES_TBI" "${LR_DIR}/${TBI_FILENAME}"

# --- Step 7: Verify all expected files ---
echo ""
echo "🧾 Verifying final file list..."
MISSING=0
while IFS= read -r fname; do
    if [ ! -f "${LR_DIR}/${fname}" ]; then
        echo "  ❌ Missing: ${LR_DIR}/${fname}"
        MISSING=$((MISSING + 1))
    fi
done <<< "$FILENAMES"

if [ "$MISSING" -eq 0 ]; then
    ok "All files accounted for."
else
    warn "$MISSING file(s) still missing after extraction."
fi

# --- Step 8: Generate .fai index for reference sequences ---
echo ""
echo "🧬 Generating .fai index files for reference sequences (Haplotype == 'ref')..."
INDEXED=0
while IFS=$'\t' read -r hap fasta; do
    FASTA_PATH="${LR_DIR}/${fasta}"
    FAI_PATH="${FASTA_PATH}.fai"
    if [ ! -f "$FASTA_PATH" ]; then
        warn "Skipping $fasta — file not found."
        continue
    fi
    if [ -f "$FAI_PATH" ]; then
        ok "$FAI_PATH already exists."
        INDEXED=$((INDEXED + 1))
        continue
    fi
    log "Indexing $fasta ..."
    if samtools faidx "$FASTA_PATH"; then
        ok "Created $FAI_PATH."
        INDEXED=$((INDEXED + 1))
    else
        warn "Failed to index $fasta."
    fi
done < <(python3 -c "
import csv, sys
with open('${CSV_PATH}', newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        hap = row.get('Haplotype','').strip().lower()
        fn  = row.get('Filename','').strip()
        if hap == 'ref' and fn:
            print(hap + '\t' + fn)
")

ok "Indexed $INDEXED reference file(s)."

echo ""
echo "🎉 Setup complete! Activate your environment with: mamba activate invntr / conda activate invtr"
echo ""
