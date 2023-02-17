#!/usr/bin/env python
#-*- coding: UTF-8 -*-
from __future__ import division
import sys
import argparse
import operator
from tabulate import tabulate
import re
import plotext as plt
program = 'fastaMetrics.py'

parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description="\n\n\tCalculate summary statistics for genome assemblies (multifasta file)\n---------------------\npython 3.0 or higher\n\n")
requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-i", "--input", dest='input', required=True, type=str, help='Input fasta file (or genome assembly)')
requiredNamed.add_argument("-o", "--output", dest='output', required=True, type=str, help='Prefix: Output text file, results are also printed on screeen')
args = parser.parse_args()

#requiredNamed = parser.add_argument_group('required named arguments')
#requiredNamed.add_argument('-i', '--input', help='Input file name', required=True)

sequence = ""
fasta_dir = {}
fas_len = {}
gc_dir = {}
lens_f = []
seq_gc = ''

def count_Ns(f_sequence):
	nseq = f_sequence.upper()
	totns = []
	single_n = 0
	for match in re.finditer('N+', str(nseq)):
		totns.append(match.end() - match.start())
	for val in totns:
		if val <= 5:
			single_n += 1
		else:
			continue
	ns_in_seq = sum(totns)
	gaps = len(totns) - single_n    #gaps = n>5
	#Return a list of values [total Ns in fasta, gaps(N strings longer to 5bp)]
	return [ns_in_seq, gaps]

total_Ns = []                                     #-------->
total_gaps =[]                                    #-------->
with open(args.input, 'r') as f:
	for line in f:
		if line.startswith('>'):
			if sequence:
				single_gc = float((sequence.count('G') + sequence.count('C')))/len(sequence) * 100
				fas_len[seq_id] = len(sequence)
				gc_dir[seq_id] = single_gc		
				fasta_dir[seq_id] = sequence
				lens_f.append(len(sequence))
				seq_gc += sequence
				n_metrics = count_Ns(sequence)     #-------->
				total_Ns.append(n_metrics[0])      #-------->
				total_gaps.append(n_metrics[1])    #-------->
				sequence = ""
			seq_id = line.rstrip()[1:]
		else:
			sequence = sequence + line.strip()
	fasta_dir[seq_id] = sequence
	lens_f.append(len(sequence))
	seq_gc += sequence
	gc_dir[seq_id] = single_gc
	fas_len[seq_id] = len(sequence)
	n_metrics = count_Ns(sequence)     #-------->
	total_Ns.append(n_metrics[0])      #-------->
	total_gaps.append(n_metrics[1])    #-------->
gc_cont = float((seq_gc.count('G') + seq_gc.count('C')))/len(seq_gc) * 100

def calculate_N50(list_of_lengths):
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return median
n50_res = calculate_N50(lens_f)
mean_len = sum(lens_f)/len(lens_f)

with open(args.output + '.txt', 'w') as out:
	out.write('\nTotal sequences in fasta:\t' + str(len(fasta_dir.keys()))+'\n')
	out.write('Total bases in fasta:\t' + str(len(seq_gc))+'\n')
	out.write('Total Ns in fasta:\t' + str(sum(total_Ns))+'\n')                                        #-------->
	out.write('Total Gaps in fasta (N strings longer to 5):\t' + str(sum(total_gaps))+'\n')            #-------->
	out.write('Shortest sequence in fasta file:\t' + str(min(fas_len.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + min(fas_len.items(), key=operator.itemgetter(1))[0]+')\n')
	out.write('Longest sequence in fasta file:\t' + str(max(fas_len.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + max(fas_len.items(), key=operator.itemgetter(1))[0]+')\n')
	out.write('Mean length:\t' + str(mean_len)+'\n')
	out.write('N50:\t' + str(n50_res)+'\n')
	out.write('Overall G+C content:\t' + str(gc_cont)+'\n')
	out.write('Min G+C per sequence in fasta file:\t' + str(min(gc_dir.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + min(gc_dir.items(), key=operator.itemgetter(1))[0]+')\n') 
	out.write('Max G+C per sequence in fasta file:\t' + str(max(gc_dir.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + max(gc_dir.items(), key=operator.itemgetter(1))[0]+')')

final_stats = {}
final_stats['a. Total sequences'] = str(len(fasta_dir.keys()))
final_stats['b. Total bases'] = str(len(seq_gc))
final_stats['c. Total Ns'] = str(sum(total_Ns))                                 #-------->
final_stats['d. Total Gaps (N strings longer to 5)'] = str(sum(total_gaps))     #-------->
final_stats['e. Shortest sequence'] = str(min(fas_len.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + min(fas_len.items(), key=operator.itemgetter(1))[0] + ')'
final_stats['f. Longest sequence'] = str(max(fas_len.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + max(fas_len.items(), key=operator.itemgetter(1))[0]+ ')'
final_stats['g. Mean length'] = str(sum(lens_f)/len(lens_f))
final_stats['h. N50'] = str(n50_res)
final_stats['i. Min G+C'] = str(min(gc_dir.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + min(gc_dir.items(), key=operator.itemgetter(1))[0]+ ')'
final_stats['j. Max G+C'] = str(max(gc_dir.items(), key=operator.itemgetter(1))[1]) + ', sequence name(' + max(gc_dir.items(), key=operator.itemgetter(1))[0]+ ')'
final_stats['k. Overall G+C content'] = str(gc_cont)
final_stats['l. Results available at'] = str('./' + args.output + '_Summary.txt, ' + args.output + '_CumulativeLength.txt')

print ('\n\n')
print ('\t\t----------------- Fasta statistics result -----------------')
print ((tabulate(sorted(final_stats.items()), headers=['Stats','Value'], tablefmt="fancy_grid")))
print ('\n\n')
print ('\n\t\t----------------- Plotext Output -----------------\n')
#lens_f
cum_list=[]
y = 0 
ct_cum = 0

len_sorted = sorted(lens_f, reverse=True)
with open(args.output + '_CumulativeLength.txt', 'w') as out_c:
	out_c.write('Contig_Number\tCumulativeLength\n')
	for x in range(0,len(len_sorted)):
		ct_cum += 1
		y+=len_sorted[x]
		cum_list.append(y)
		out_c.write(str(ct_cum) + '\t' + str(y) + '\n')
plt.plot(cum_list)
plt.scatter(cum_list)
plt.figsize(60, 20)     #######para la version mas reciente de python "plotsize" debe reemplazarse por "plot_size"
plt.grid(True)
plt.title("Cumulative length of fasta sequences")
plt.xlabel("Contig_number")
plt.show()





