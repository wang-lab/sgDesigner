# sgDesigner
WashU gRNA Designer for CRISPR/Cas9 Knockout

Welcome to sgDesigner, the Washington University gRNA designer for CRISPR/Cas9 knockouts. This program is distributed under the GNU General Public License as published by the Free Software Foundation. You may use it freely in your research with appropriate acknowledgment. The program is provided "as is" and the authors are not responsible for consquences from the use of the program.

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
* The versions of Python packages which sgDesigner used are, specifically:

```
   NumPy: 1.15.2
   SciPy: 1.1.0
   Scikit-learn: 0.20.0
   XGBoost: 0.80
```
   * [Installation instruction of Python packages](https://packaging.python.org/tutorials/installing-packages/)
   
## Installation of sgDesigner standalone program

1. Place the sgDesigner.tar.gz file anywhere in your Linux system and uncompress using the following command:

```
   'tar -xzvf sgDesigner.tar.gz'
```
   The installation should complete within 5 seconds.

2. Copy your input FASTA files into the newly created sgDesigner directory.

3. Type 'perl sgDesigner.pl' to run the program and view the help file.

A README file is included in the sgDesigner standalone package, with examples and detailed explanations of the commands available.

