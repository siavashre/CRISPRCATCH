import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input Cycle file", required=True)
parser.add_argument("-o", "--output", help="Output directory", required=True)
args = parser.parse_args()
cycle_dir = args.input
output_dir = args.output

def print_cycle(cycle_path, cycle_num, seg_dict):
	segment_counter = {}
	segment_counter = defaultdict(lambda: 0, segment_counter)
	for i in cycle_path.strip().split(','):
		seg_id = int(i[:-1])
		if seg_id!=0:
			segment_counter[seg_id]+=1
	with open(output_dir+'_cycle'+cycle_num+'.bed', 'w') as f :
		for k in seg_dict.keys():
			if segment_counter[k] > 0:
				line = "{chrom}\t{start}\t{end}\t{cp}\n".format(chrom = seg_dict[k][0], start = seg_dict[k][1], end = seg_dict[k][2], cp = segment_counter[k])
				f.write(line)


seg_dict = {}
with open(cycle_dir,'r') as f:
	for line in f :
		if line.startswith('Segment'):
			line = line.strip().split('\t')
			seg_id = int(line[1])
			seg_chr = line[2]
			seg_start = int(line[3])
			seg_end = int(line[4])
			seg_dict[seg_id] = [seg_chr, seg_start, seg_end]
		elif line.startswith('Cycle'):
			line = line.strip().split(';')
			cycle_num = line[0].split('=')[1]
			cycle_path = line[3].split('=')[1]
			print_cycle(cycle_path, cycle_num,seg_dict)
