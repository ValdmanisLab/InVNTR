# InVNTR
Isolate and neaten VNTRs (InVNTR): A command-line tool for VNTR analysis and visualization.
### Details:
InVNTR allows the rapid extraction of tandem repeats from consensus sequences, like genomic fasta. We use it with high quality long read assemblies, such as those available from the HPRC, HGSVC, and the 1000 Genomes Project. Rather than a reference genome, we use something like in-silico PCR in order to isolate the loci of interest. For this reason, VNTR need to have conserved and unique flanking sequence on at least one side. InVNTR assumes all user arguments have the same strandedness, but will search for reverse complements in fasta if it cannot find sequences that match the arguments.

## Usage: To isolate VNTR alleles and neaten motifs from genomes.
### Input: 
*Required:*
- --start | conserved and unique starting sequence, either first motif of the VNTR or sequence directly before the VNTR.
- --end **OR** --no_end | if you use the end flag, it is used in the same way as the start flag, by providing a conserved and unique ending sequence, either the last motif of the VNTR or the sequence directly following the end of the VNTR. If you use the no_end flag, then you need to specify a numerical value for the sequence to put the end for. This can be useful if the sequence on one end of the VNTR is not conserved.
- --folder | directory containing .fasta or .fa files
- --length **OR** --delimiter | if your VNTR has a consistent motif length, you may provide --length as a numerical value for the number of bases in a single repeat motif. If your VNTR has variable length you can provide a delimiter, a set of characters that appear at the beginning of each repeat motif.

*Optional:*
- --after and --before | These variables(y or n, default n) let you specify if the sequence that you used for start and end are before or after the actual VNTR, in case the start of the VNTR is variable or repeated in the VNTR. If you specify y, then InVNTR will remove the length of the start or end from the VNTR in the outputs.
- --type | filetype (default is .fa or .fasta), which allows you to operate InVNTR on directories which contain other kinds of files.
- --name | (name to append to output files, default is VNTR)
- --max_allele_length | Maximum nucleotide length that will be written to the alleles file output, per file. Default is 10,000. This stops the alleles file from getting too large if there is an issue with the start and end.


## Dependencies:
### Python 3
*non-built-in modules:*
- Bio pandas argparse
- Which can be installed with pip on Windows, Linux, or Mac:
```
 python3 -m pip install Bio pandas argparse 
```

### For MacOS

For MacOS you may need the homebrew package c-blosc and the python module tables before you can install the above required packages.

```
brew install c-blosc
python3 -m pip install tables
```

## Examples
Outputs are saved in the directory the command is run from. Execute the command in the terminal from the directory you wish them to be saved in. Outputs will be overwritten if given the same name as a previous output.
### Windows Motif Length Example: 
~~~~
python C:\path_to\InVNTR_1.4.py -f E:\wholegenomes -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL -l 30
~~~~
<details>
<summary>More detailed explanation</summary>
This example extracts the IVL VNTR. The IVL VNTR is ideal for the motif length option because it has a consistent motif length of 30. In this example the start is the first two motifs, as together they are consistent accross alleles, but the sequence in unique in the genome. The end is the last motif.  

</details> 
<br />

### Windows Delimiter Example: 
~~~~
python C:\path_to\InVNTR_1.4.py -f E:\wholegenomes -s CATCTCCTCCTCCTCACCTCCTGCTGTGGTGCACAGATACCTATAGGCAGGCTC -e CATCTCCTCCTCCCGAGCTCCTCCCCTAGTGCACAGATACCTATAGGCAGGCTC -n SORL1 -d CATCT
~~~~
<details>
<summary>More detailed explanation</summary>
This example extracts the SORL1 VNTR. The SORL1 VNTR is ideal for the delimiter option because it has a consistent sequence at the start of each motif, despite variable length.

</details> 
<br />

### Unix (Mac/Linux) Example: 
~~~~
python ~/path_to/InVNTR_1.4.py -f /Volumes/wholegenomes -l 30 -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL
~~~~
## Outputs:
- VNTR.csv – excel file with alleles on the x axis, with each allele named after the file and the tandem repeat cut into motifs based on the length or delimiter provided from the beginning of the allele.
<br />
- VNTR_allele_length.csv – excel file with alleles on the y axis, with each allele named after the file length of the VNTR for each allele.
<br />
- VNTR_alleles.txt – text file with each allele after the name of the file.
<br />
- VNTR_errors.txt – text file with any errors (ie. '[filename] has only reverse end’ or ‘[filename] has no sign of tandem repeat in forward, reverse’).
<br />
- VNTR_frequency.csv – csv with motifs sorted by frequency accross the entire population. 
<br />

## Current help flag output:

usage: InVNTR_1.4.py [-h] -f  [-l ] [-d ] -s  [-e ] [-ne ] [-mal ] [-b ] [-a ] [-t ] [-n ]

Extract VNTR

options:
  -h, --help            show this help message and exit <br />
  -f , --folder         Path to folder. ex: E:\HPRC <br />
  -l [], --length []    length of consensus motif. ex: 30 : DO NOT USE IF USING DELIMITER, LENGTH TAKES PRIORITY <br />
  -d [], --delimiter []
                        consistent beginning of motif: DO NOT USE IF USING LENGTH, LENGTH TAKES PRIORITY <br />
  -s , --start          beginning of VNTR. ex: GCTA <br />
  -e [], --end []       end of VNTR. ex: GCTA <br />
  -ne [], --no_end []   If there is no end, when should the program cut out? default is 10000 <br />
  -mal [], --max_allele_length []
                        Maximum nucleotide length that will be written to the allele file output, per file. Default is 10000 <br />
  -b [], --before []    Is the start value of the VNTR before the VNTR? y OR n <br />
  -a [], --after []     Is the end value of the VNTR after the VNTR? y OR n <br />
  -t [], --type []      Filetype ex: .fasta (default is .fa or fasta) <br />
  -n [], --name []      output name) <br />

## Future Direction:
Improving on motif decomposition remains the top priority. Refactoring is on the horizon as well.