# sgDesigner
WashU gRNA Designer for CRISPR/Cas9 Knockout

Welcome to sgDesigner, the Washington University gRNA designer for CRISPR/Cas9 knockouts. This program is distributed under the GNU General Public License as published by the Free Software Foundation. You may use it freely in your research with appropriate acknowledgment. The program is provided "as is" and the authors are not responsible for consquences from the use of the program.

A README file is included in the sgDesigner standalone package, with examples and detailed explanations of the commands available.

## Requirements

* This package is supported for *Linux* operating systems. The package has been tested on the following systems:
```
   Linux: CentOS 7.1.1503
```
* Perl 5 interpreter or higher on a Red-Hat compatible Linux system is required.
   * [Installation instruction](https://learn.perl.org/installing/)
   * [Perl download](https://www.perl.org/get.html)
* Python 3.7.0 is used to generate the model and do the prediction. 
   * [Installation instruction](https://realpython.com/installing-python/)
   * [Python download](https://www.python.org/downloads/)
   * [Installation instruction of Python packages](https://packaging.python.org/tutorials/installing-packages/)
* The versions of Python packages which sgDesigner used are, specifically:
```
   NumPy: 1.15.2
   SciPy: 1.1.0
   Scikit-learn: 0.20.0
   XGBoost: 0.80
```
  
   
## Installation of sgDesigner standalone program

* Place the sgDesigner.tar.gz file anywhere in your Linux system and uncompress using the following command:
```
   tar -xzvf sgDesigner.tar.gz
```
* Copy your input FASTA files into the newly created sgDesigner directory.
* Type 'perl sgDesigner.pl' to run the program and view the help file.


## Command Line Parameters

* Direct sequence submission (-s or --sequence):
   
   This option allows the user to submit a single sequence directly for analysis, using the following command:
   ```
      perl sgDesigner.pl –s <sequence>
   ```
   This option is most useful for users who wish to determine the efficacy of a single gRNA. Any sequences submitted must be at least 26 bases long (including the NGG PAM region) and contain only A, T, U, C, or G. These rules also apply for any FASTA sequences that are submitted, which are covered in more detail in the next section.
* FASTA file submission (-f or --file):
   
   This option allows the user to submit one or more sequences in a FASTA file, using the following command:
   ```
      perl sgDesigner.pl –f myFastaFile.fasta
   ```
   This should be provided in FASTA format. In a FASTA file, a definition line that begins with begins with ‘>’ is required for each DNA sequence. For example:
   ```
      >gi|4507798|ref|NM_000462.1| Homo sapiens ubiquitin protein ligase E3A (UBE3A), mRNA
      ATGGAGAAGCTGCACCAGTGTTATTGGAAATCAGGAGAACCTCAGTCTGACGACATTGAAGCTAGCCGA
      TGAAGCGAGCAGCTGCAAAGCATCTAATAGAACGCTACTACCACCAGTTAACTGAGGGCTGTGGAAATA
      AGCCTGCACGAATGAGTTTTGTGCTTCCTGTCCAACTTTTCTTCGTATGGATAATAATGCAGCAGCTAT
      TAAAGCCCTCGAGCTTTATAAGATTAATGCAAAACTCTGTGATCCTCATCCCTCCAAGAAAGGAGCAAG
      CGCAGCTTACCTTGAGAACTCGAAAGGTGCCCCCAACAACTCCTGCTCTGAGATAAAAATGAACAAGAA
      AGG
   ```
   Submitted sequences must be between 26 and 100,000 nt in length and contain A, T, U, C, or G. Three sample files are also provided: one containing a single short sequence (30 nt), one with a single long sequence (8,322 nt), and one with 3 sequences of 300, 600, and 300 bases long.
* Sample file submission:
   This option allows the user to try one of the three previously mentioned sample files, using the following command:
   ```
      perl sgDesigner.pl –e <short|long|multiple>
   ```
   One of the three options shown will call the respective sample file.
   
## Sample expected outputs
