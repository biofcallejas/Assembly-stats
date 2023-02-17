# Genome Assembly Stats
## Python script for summarizing a genome assembly

This script will summaryze the nucleotide content in a fasta file by calculating:

```
Total sequences
Total bases
Total Ns (undefined bases)
Total Gaps (strings of undefined bases longer than 5)
Shortest sequence
Longest sequence
Mean length
N50
Overall G+C content
Min G+C per sequence
Max G+C per sequence
```
## Running the script with test files

Before running the script, make sure the following libraries are installed (tested in Python 3.8.12):

```
import tabulate
import plotext as plt
```
Results are also printed on the screen, make sure you're using a suitable python version.

```
python fastaMetrics_py3.py -h 
usage: fastaMetrics.py [-h] -i INPUT -o OUTPUT

Calculate summary statistics for genome assemblies
---------------------
Tested in python 3.8.12

optional arguments:
  -h, --help        show this help message and exit

Mandatory arguments:
  -i INPUT, --input, Input fasta file (genome assembly)
  -o OUTPUT, --output, Prefix for the output files, results are also printed on screeen
```
Running script with test file:

```
python fastaMetrics_py3.py -i pfalciparium3D7_PlasmoDB-61.fasta -o pf3d7
```

![image](https://user-images.githubusercontent.com/45425927/219796445-5bcce619-68b8-44d5-927a-ac0391bf8917.png)


## Output

*In process*
