from collections import defaultdict
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input bed file", required=True)
parser.add_argument("-o", "--output", help="Input bed file", required=True)
args = parser.parse_args()
l = []
with open(args.input, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        cnv = int(line[3])
        l.append(cnv)

mean = int(np.mean(l) * 2.5)
# with open(args.input[:-4]+'_norm.bed', 'w') as g:
with open(args.output, 'w') as g:
	with open(args.input, 'r') as f:
	    for line in f:
	        line = line.strip().split('\t')	      
	        cnv = int(line[3])
	        line[3] = str(min(cnv,mean))
	        line2 = '\t'.join(i for i in line) + '\n'
	        g.write(line2)
