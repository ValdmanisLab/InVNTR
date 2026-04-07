# InVNTR
Isolate and Neaten VNTRs (InVNTR): A command-line tool for VNTR analysis and visualization.

### Details:
InVNTR allows the rapid extraction of tandem repeats from consensus sequences, like genomic FASTA. We use it with high quality long read assemblies, such as those available from the HPRC, HGSVC, and the 1000 Genomes Project. Rather than a reference genome, we use something like in-silico PCR in order to isolate the loci of interest. For this reason, VNTRs need to have conserved and unique flanking sequence. 

---

## What's New in Version 1.5 

- **Coordinate mode (now default):** You can now specify a VNTR locus using genomic coordinates in interval format (e.g. `chr12:40,479,053-40,491,305` for the VNTR exon of MUC19) instead of manually providing `--start` and `--end` sequences. InVNTR will automatically extract flanking sequences from a reference genome. GRCh38 is the default reference, though T2T coordinates should also work. InVNTR relies on .vcf files to generate wildcards in flanking regions for snp's in the autosome and X chromosome. InVNTR is currently lacking a y chromosome .vcf, so y chromosome coordinates are not currently usable.
- **`--start` and `--end` are now optional** when using coordinate mode. They remain available if you prefer to provide sequences manually.
- **Assemblies file is now required** for coordinate mode. See the Setup section below.
- **Default output is now multi-sheet excel workbook.** See the updated Output section below.

---

## Setup

Before running InVNTR 1.5, you need to download the required reference assemblies. A setup script is provided to automate this, as well as configure an environment for running InVNTR. You must already have mamba or conda installed and active on your terminal application, as well as wget.

### 1. Download the setup script

```
wget https://raw.githubusercontent.com/ValdmanisLab/InVNTR/main/Stable/setup_InVNTR_1.5.sh
```

### 2. Run the setup script

Run from the directory where you want the assemblies to be downloaded:

Full dataset - 4.25 terabytes, expect setup to take no less than 3 days.
```bash
# Download full assembly manifest and all files
./setup_InVNTR_1.5.sh
```

Testing dataset - ~100GB terabytes, expect setup to take no less than 1 hour.
```bash
./setup_InVNTR_1.5.sh -t
# OR
./setup_InVNTR_1.5.sh -testing
```

The `-t` / `--testing` flag uses `assemblies_testing.csv`, which contains a reduced set of samples — useful for testing your setup before committing to the full download.

The setup script will:
1. Download the appropriate assembly manifest CSV from the InVNTR GitHub repo.
2. Setup an 'invntr' conda environment with all of the dependencies for decompressing the assemblies and running invntr.
3. Check for missing files and download our .agc from our sharepoint and extract them.
5. Generate `.fai` index files for reference sequences using `samtools faidx`.
6. Download a 1000 Genomes Project phased VCF and it's index.

---

## Usage: To Isolate VNTR Alleles and Neaten Motifs from Genomes

### Input:

*Required:*
- `--folder` | directory containing `.fasta` or `.fa` files
- `--length` **OR** `--delimiter` | if your VNTR has a consistent motif length, provide `--length` as a number of bases per motif. If your VNTR has variable-length motifs, provide `--delimiter`, a set of characters that appear at the start of each motif.
- `--coordinate` **OR** `--start` | in coordinate mode (default), provide a genomic interval like `chr12:40,479,053-40,491,305`. Alternatively, use `--start` to provide the conserved sequence at the beginning of the VNTR directly.

*When using `--start` manually:*
- `--end` **OR** `--no_end` | provide a conserved ending sequence, or use `--no_end` with a numerical cutoff if the end is not conserved.

*Optional:*
- `--reference` / `-r` | reference genome to use with coordinate mode (default: GRCh38; T2T should also work)
- `--flank` | number of bp for flanking sequences extracted around the coordinate interval (default: 50)
- `--after` and `--before` | (`y` or `n`, default `n`) specify whether your `--start` or `--end` sequences are outside the VNTR. If `y`, InVNTR removes the length of those sequences from the output.
- `--type` | filetype to look for in the folder (default: `.fa` or `.fasta`)
- `--name` | prefix for output files (default: `VNTR`)
- `--max_allele_length` | maximum nucleotide length written to the alleles file per sample (default: 10,000)
- `--print` / `-p` | if set, exports each sheet of the Excel output as individual CSV or TXT files
- `--verbose` / `-v` | enables verbose output

---

## Examples

Outputs are saved in the directory the command is run from. Execute in the terminal from the directory you want outputs saved to. InVNTR will append a _# rather than overwrite a previous output.

### Coordinate Mode (New Default) — Unix:
```bash
python3 ~/path_to/InVNTR_1.5.py -f /Volumes/wholegenomes -c chr1:152910255-152911421 -l 30 -n IVL
```
<details>
<summary>More detailed explanation</summary>
This uses coordinate mode to extract the IVL VNTR from GRCh38 coordinates. InVNTR automatically derives flanking sequences from the reference, so no --start or --end is needed. The motif length of 30 is still required to decompose the repeat.
</details>
<br />

### Coordinate Mode — Windows:
```
python C:\path_to\InVNTR_1.5.py -f E:\wholegenomes -c chr1:152910255-152911421 -l 30 -n IVL
```
<br />

### Manual Start/End Mode — Windows Motif Length Example:
```
python C:\path_to\InVNTR_1.5.py -f E:\wholegenomes -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL -l 30
```
<details>
<summary>More detailed explanation</summary>
This example extracts the IVL VNTR using manually provided start and end sequences. The IVL VNTR is ideal for the motif length option because it has a consistent motif length of 30. The start is the first two motifs, which together are consistent across alleles and unique in the genome. The end is the last motif.
</details>
<br />

### Manual Start/End Mode — Windows Delimiter Example:
```
python C:\path_to\InVNTR_1.5.py -f E:\wholegenomes -s CATCTCCTCCTCCTCACCTCCTGCTGTGGTGCACAGATACCTATAGGCAGGCTC -e CATCTCCTCCTCCCGAGCTCCTCCCCTAGTGCACAGATACCTATAGGCAGGCTC -n SORL1 -d CATCT
```
<details>
<summary>More detailed explanation</summary>
This example extracts the SORL1 VNTR using a delimiter. The SORL1 VNTR is ideal for the delimiter option because it has a consistent sequence at the start of each motif, despite variable motif length.
</details>
<br />

### Manual Start/End Mode — Unix (Mac/Linux) Example:
```bash
python3 ~/path_to/InVNTR_1.5.py -f /Volumes/wholegenomes -l 30 -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL
```

---

## Outputs

The primary output is a single `.xlsx` Excel workbook (default: `VNTR.xlsx`). If a file with that name already exists, a counter is appended (e.g. `VNTR_1.xlsx`). It contains the following sheets:

**SeattlePlot** — The main data matrix. Columns are individual genome assemblies (named after their source file), rows are motif positions (1, 2, 3...). Each cell contains the sequence of that motif at that position in that assembly. The first four rows below the header contain per-sample metadata: Allele Length, Population code, Superpopulation code, and File Origin.

**Allele_Length** — One row per assembly with columns: Filename, Allele Length (bp), Population code, Superpopulation code, and File Origin.

**Frequency** — One row per unique motif sequence observed across all samples and positions, with columns: Motif, Frequency (count), and Length (bp). Useful for identifying dominant vs. rare motif variants.

**Alleles** — FASTA-formatted records of each assembly's full extracted VNTR sequence (truncated to `--max_allele_length` if needed), stored as alternating `>filename` / sequence rows.

**Population** — Statistical summary of allele length variation. Contains mean and standard deviation grouped by Population code, then by Superpopulation code. Also includes a one-way ANOVA testing whether allele lengths differ significantly between groups, followed by a Tukey HSD post-hoc test identifying which specific population pairs differ. Results are reported for both population-level and superpopulation-level groupings.

**Errors** — Two columns (Sample, Error) listing any assemblies that failed to extract, with the reason (e.g. missing flanking sequences, partial matches, file read errors).

**Log** — Timestamped record of every processing step, the full command used to invoke the script, total elapsed time, number of genomes processed, and average time per genome.

---

If `--print` / `-p` is set, each sheet is additionally exported as a separate file alongside the `.xlsx`: SeattlePlot, Allele_Length, Frequency, and Population as `.csv`; Alleles, Errors, and Log as `.txt`. All files are prefixed with the `--name` argument.

---

## Future Direction

Combining .vcf files and adding a y chromosome .vcf are the current top priority. Adding Alu and other RE detection and avoidance are also a priority. As well as continuing to think about repeat motif decomposition. The SeattlePlot function is currently underutilized due to his rigidity.
