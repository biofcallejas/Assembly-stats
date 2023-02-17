# Genome Assembly Stats
## Python script for summarizing a genome assembly

Before running the script, make sure the following libraries are installed (tested in Python 3.8.12):

```
import tabulate
import plotext as plt
```
## Running the script with test files

```
python fastaMetrics_py3.py -h 
usage: fastaMetrics.py [-h] -i INPUT -o OUTPUT

Calculate summary statistics for genome assemblies (multifasta file)
---------------------
Tested in python 3.8.12

optional arguments:
  -h, --help        show this help message and exit

Mandatory arguments:
  -i INPUT, --input, Input fasta file (genome assembly)
  -o OUTPUT, --output, Prefix for the output files, results are also printed on screeen
```
## Output

*In process*
