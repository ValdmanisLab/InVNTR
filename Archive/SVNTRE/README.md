# InVNTR
Isolate and neaten VNTRs (InVNTR): A command-line tool for VNTR analysis and visualization.

# Updated documentation for 1.4 incoming

## Usage: To isolate VNTR alleles and neaten motifs from genomes.
### Input: 
*Required:*
- consistent and unique starting and ending sequence (either first & last motifs or directly outside of them)
- folder containing .fasta or .fa files (whole genome or part)
- motif length **OR** delimiter - (recommend motif length if your motifs are of a consistent length, delimiter if they vary in length but all begin with the same unique sequence)

*Optional:*
- after and before variable (program defaults to assuming you are using the first and last motif, if, instead you are using sequence directly outside the VNTR, you must specify so that it is not included in the extraction)
- filetype (.fa or .fasta only, default is .fa)
- name (name to append to output files, default is VNTR)


## Dependencies:
### Python 3
*non-built-in modules:*
- Bio pandas argparse
- Which can be installed with pip:\
python3 -m pip install Bio pandas argparse

### MacOS

Additionally, for MacOS you will need the module tables, and the homebrew package c-blosc before you can install the above required packages. Probably.

## Examples
### Windows Motif Length Example: 
Execute this in terminal from the directory you want the output to be in:
~~~~
C:/Users/username/AppData/Local/Programs/Python/Python310/python.exe C:/VNTR/InVNTR_1.4.py -f E:\wholegenomes -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL -l 30
~~~~
<details>
<summary>More detailed explanation</summary>
This example extracts the IVL VNTR. The IVL VNTR is ideal for the motif length option because it has a consistent motif length of 30. In this example the start is the first two motifs, as together they are consistent accross alleles, but the sequence in unique in the genome. The end is the last motif.  

</details> 
<br />

### Windows Delimiter Example: 
Execute this in terminal from the directory you want the output to be in:
~~~~
C:/Users/username/AppData/Local/Programs/Python/Python310/python.exe C:/VNTR/InVNTR_1.4.py -f E:\wholegenomes -s CATCTCCTCCTCCTCACCTCCTGCTGTGGTGCACAGATACCTATAGGCAGGCTC -e CATCTCCTCCTCCCGAGCTCCTCCCCTAGTGCACAGATACCTATAGGCAGGCTC -n SORL1 -d CATCT
~~~~
<details>
<summary>More detailed explanation</summary>
This example extracts the SORL1 VNTR. The SORL1 VNTR is ideal for the delimiter option because it has a consistent sequence at the start of each motif, regardless of its length.

</details> 
<br />

### Mac Example: 
Execute this in terminal from the directory you want the output to be in which also contains the .py script file:
~~~~
Python3 InVNTR_1.4.py -f /Volumes/wholegenomes -l 30 -s TGGAGCTCCCAGAGCAGCAGGAGGGGCACCTGAAGCACCTAGAGCAGCAGGAGGGACAGC -e TGGAGCAGCAGAAGGGGCAGCTGGAGCAGC -n IVL
~~~~
## Outputs:
- VNTR.csv                       – excel file with each allele named after the file, cut into simple motifs based on the length or delimiter provided from the beginning of the allele
<br />

- VNTR_alleles.txt                   – text file with each allele after the name of the file.
<br />
- VNTR_errors_and_frequency.txt      – text file with any errors (ie. ‘Filename has only reverse end’ or ‘filename has no sign of tandem repeat in forward, reverse’) and ranked frequency of simple motifs.


## Current help flag output:

usage: InVNTR_1.4.py [-h] -f  [-l ] [-d ] -s  -e  [-b ] [-a ] [-t ] [-n ]

Extract VNTR

options:
  -h, --help            show this help message and exit<br />
  -f , --folder         Path to folder. ex: E:\HPRC<br />
  -l [], --length []    length of consensus motif. ex: 30 : DO NOT USE IF USING DELIMITER, LENGTH TAKES PRIORITY<br />
  -d [], --delimiter [] consistent beginning of motif: DO NOT USE IF USING LENGTH, LENGTH TAKES PRIORITY<br />
  -s , --start          beginning of VNTR. ex: GCTA<br />
  -e , --end            end of VNTR. ex: GCTA<br />
  -b [], --before []    Is the start value of the VNTR before the VNTR? y OR n<br />
  -a [], --after []     Is the end value of the VNTR after the VNTR? y OR n<br />
  -t [], --type []      Filetype ex: .fasta (as default is .fa)<br />
  -n [], --name []      output name)<br />

## Future Direction:
Improving on motif decomposition remains the top priority.