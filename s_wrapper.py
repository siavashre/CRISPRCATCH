import os
from collections import defaultdict
import pandas as pd
import argparse
import csv
#####	Helper Functions ###########
class Discordant_edge:
	chr1 = ''
	chr2 = ''
	pos1 = 0
	pos2 = 0
	dir1 = ''
	dir2 = ''
	cp = 0
	line = ''
	count = 0
	in_bulk = 0

def parse_csv(file):
	a = pd.read_csv(file)
	bands = list(a['band'])
	read1 = list(a['read1'])
	read2 = list(a['read2'])
	max_len = list(a['band_max_len'])
	min_len = list(a['band_min_len'])
	guide_start = list(a['guide_start'])
	guide_end = list(a['guide_end'])
	guide_chr = list(a['guide_chr'])
	d_max = {}
	d_min = {}
	d_guids = {}
	d_read = {}
	for i in range(len(bands)):
		d_max[bands[i]] = max_len[i]
		if pd.isna(max_len[i]):
			d_max[bands[i]] = min_len[i]
		d_min[bands[i]] = min_len[i]
		if pd.isna(min_len[i]):
			d_min[bands[i]] = max_len[i]
		d_guids[bands[i]] = (guide_chr[i], guide_start[i], guide_end[i])
		d_read[bands[i]] = [read1[i], read2[i]]
	return bands, d_read,d_max, d_min, d_guids

def parsing_AA_graph(graph_dir):
	l = []
	with open(graph_dir, 'r') as f:
		for line in f:
			if line.startswith('discordant'):
				line2 = line
				line = line.strip().split('\t')
				d = Discordant_edge()
				d.chr1 = line[1].split(':')[0]
				d.chr2 = line[1].split('>')[1].split(':')[0]
				d.pos1 = float(line[1].split('->')[0].split(':')[1][:-1])
				d.pos2 = float(line[1].split('->')[1].split(':')[1][:-1])
				d.dir1 = line[1].split('->')[0].split(':')[1][-1]
				d.dir2 = line[1].split('->')[1].split(':')[1][-1]
				d.cp = float(line[2])
				d.line = line2
				l.append(d)
	return l 

def compare_discordants(a,b):
	if a.chr1 == b.chr1 and a.chr2 == b.chr2 and a.dir2 == b.dir2 and a.dir1 == b.dir1 and abs(a.pos1 - b.pos1) < T and abs(a.pos2 - b.pos2) < T:
		return True
	return False

def compare_bulk_band(bulk , band):
	bulk_count = 0
	band_count = 0
	for i in range(len(band)):
		for j in range(len(bulk)):
			if compare_discordants(band[i], bulk[j]):
				band_count +=1
				band[i].in_bulk = 1
	return band , band_count
def write_yml(save_dir,bed_dir):
	with open(save_dir+'.yaml', 'w') as f:
		f.write('primary_feature_bedgraph: '+bed_dir+'\n')
		f.write('rescale_by_count: True\n')
		f.write('primary_smoothing: 500\n')
		f.write('end_trim: 150\n')
		f.write('primary_kwargs:\n')
		f.write(' linewidth: 0.1')

def contain_discordant(file):
	with open(file,'r') as f:
		for line in f:
			if line.startswith('discordant'):
				return True
	return False

def parse_cycle(file_dir,band , amplicon_number,band_max_length,band_min_length):
	rows = []
	percent_breakpoints_matched = float(match_count[band+'_amplicon'+amplicon_number]) /len(d[band+'_amplicon'+amplicon_number])
	segments = {}
	with open(file_dir,'r') as f:
		for line in f:
			if line.startswith('Segment'):
				line = line.strip().split('\t')
				seg_id = int(line[1])
				seg_chr = line[2]
				seg_start = int(line[3])
				seg_end = int(line[4])
				segments[seg_id] = (seg_chr, seg_start, seg_end)
			elif line.startswith('Cycle'):
				line = line.strip().split(';')
				in_cut_site = False
				cycle_number = line[0].split('=')[1]
				rec_path = line[3].split('=')[1].split(',')
				for seg in rec_path:
					seg = int(seg[:-1])
					if seg != 0:
						if segments[seg][0] == guids[band][0]:
							if segments[seg][1]-21 <guids[band][1]< segments[seg][2]+21:
								if segments[seg][1]-21 <guids[band][2]< segments[seg][2]+21:
									in_cut_site= True
									break
				if line[3].split('=')[1].startswith('0'):
					cyclic = 'False'
				else:
					cyclic = 'True'
				reconstruct_length = line[4].split('=')[1]
				RMSR = line[5].split('=')[1]
				DBI = line[6].split('=')[1]
				Filter = line[7].split('=')[1]
				rows.append([band, cycle_number,str(int(band_min_length))+'Kbp',str(int(band_max_length))+'Kbp', amplicon_number,reconstruct_length,percent_breakpoints_matched,DBI,RMSR,Filter,cyclic,in_cut_site, band_insert_size_mean[band], band_insert_size_std[band]])
	return rows			

def parse_insert_size(file_dir):
	with open(file_dir,'r') as f:
		lines = f.readlines()
		for i in range(len(lines)):
			line = lines[i]
			if line.startswith('## METRICS'):
				line = lines[i+2]
				line = line.strip().split('\t')
				inset_mean = float(line[5]) 
				inset_std = float(line[6])
				return inset_mean, inset_std
	return 0 , 0

def quality_report():
	header = ['band','cycle_number', 'band_min_length','band_max_length', 'AA_amplicon_id', 'reconstruction_length', 'percent_breakpoints_matched', 'DBI', 'RMSR', 'FILTER','has_cycle', 'includes cut site', 'insert_size mean', 'insert_size stdev']
	with open('report.csv', 'w') as csv_file:
		csvwriter = csv.writer(csv_file)
		csvwriter.writerow(header)
		files = os.listdir('candidate_cycles/')
		lines = []
		for f in files:
			if f.endswith('candidate_cycles.txt'):
				file_dir = 'candidate_cycles/'+f
				band = f.split('_')[1]
				amplicon_number = f.split('_')[2][8:]
				band_max_length = band_size_max[band]
				band_min_length = band_size_min[band]
				lines.extend(parse_cycle(file_dir, band, amplicon_number,band_max_length,band_min_length))
		lines = sorted(lines, key = lambda x:(x[0],x[1]) )
		csvwriter.writerows(lines)

def generated_amplicon_bed_file():
	generate_bed_cmd = 'python3 {grah_to_bed} -g {graph} --unmerged'.format(grah_to_bed ='$PreAA/scripts/graph_to_bed.py', graph = args.bulk )
	print(generate_bed_cmd)
	os.system(generate_bed_cmd)

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False	

#################################################
parser = argparse.ArgumentParser()
parser.add_argument("-sname", "--sname", help="Cell Line name", required=True)
parser.add_argument("-o", "--output", help="Output dir for saving results", required=True)
parser.add_argument("-t", "--t", help="Number of thread", required=True)
parser.add_argument("-r", "--ref", help="Reference genome version", required=True)
parser.add_argument("-bed", "--bed", help="bed file describing amplicon region", required=False)
parser.add_argument("-bulk", "--bulk", help="AA breakpoint graph file for bulk cell line", required=True)
parser.add_argument("-csv", "--csv", help="CSV input file containing band information", required=True)
parser.add_argument("-sdv", "--sdv", help="insert_sdev for running AA. Default is 8.5", required=False)
parser.add_argument("-min_sup", "--min_sup", help="Minimum sup pair reads for calling discordant edges. Default is 2", required=False)
args = parser.parse_args()

#################################################
cell_line = args.sname
output_ans = args.output
if not args.bed:
	generated_amplicon_bed_file()
	amplicon_bed_file = args.bulk[args.bulk.rfind('/')+1:args.bulk.rfind('.')]+'.bed'
else:
	amplicon_bed_file = args.bed
ref_v = args.ref
bulk_AA_graph_file = args.bulk

thread = args.t
T = 101
min_pair_support = 2
insert_sdevs = 8.5
if isfloat(args.sdv):
    insert_sdevs = float(args.sdv)
if isfloat(args.min_sup):
    min_pair_support = int(args.min_sup)
########################################
band_list , band_read, band_size_max, band_size_min, guids = parse_csv(args.csv)
amplicon_mapping = defaultdict(list)
os.chdir(output_ans)
new_dir_cmd = 'mkdir graph_files'
os.system(new_dir_cmd)
fil_dir_cmd = 'mkdir filtered_graph_files'
os.system(fil_dir_cmd)
band_insert_size_mean = {}
band_insert_size_std = {}

for band in band_list:
	name = cell_line + '_' + band
	f1 = band_read[band][0]
	f2 = band_read[band][1]
	Pre_AA_cmd = "python3 {PrepareAA} --ref {ref_v} -t {thread} -s {name} --fastqs {f1} {f2} --no_filter --cngain 0 --cnsize_min 0 --cnv_bed {amplicon_bed_file} ".format(
	PrepareAA = '$PreAA/PrepareAA.py' , ref_v= ref_v, thread = thread, name = name , f1 = f1 , f2 = f2, amplicon_bed_file = amplicon_bed_file)
	print(Pre_AA_cmd)
	os.system(Pre_AA_cmd)
	os.mkdir(name+'_AA_results')
	amplified_amplicon_cmd = '$AA_SRC/amplified_intervals.py --no_cstats --bed {bed} --bam {bam} --ref hg19 --cnsize_min 0 --gain 0 --out {out}'.format(bed=amplicon_bed_file,bam = name + '.cs.rmdup.bam', out= name + '_AA_CNV_SEEDS' )
	print(amplified_amplicon_cmd)
	os.system(amplified_amplicon_cmd)
	#######		Insert Size Distribution 	########
	bam_file = name+'.cs.rmdup.bam'
	insert_size_txt = name + '_AA_results/insert_size.txt'
	insert_size_pdf= name + '_AA_results/insert_size.pdf'
	insert_size_cmd = "java -jar $PICARD/picard.jar CollectInsertSizeMetrics I={bam_file} O={insert_size_txt} H={insert_size_pdf} M=0.5".format(bam_file = bam_file , insert_size_txt = insert_size_txt , insert_size_pdf= insert_size_pdf )
	print(insert_size_cmd)
	os.system(insert_size_cmd)
	i_mean, i_std = parse_insert_size(insert_size_txt)
	band_insert_size_mean[band] = i_mean
	band_insert_size_std[band] = i_std
	AA_cmd = "$AA_SRC/AmpliconArchitect.py --out {out} --downsample -1 --bed {bed} --bam {bam} --ref hg19 --pair_support_min {min_sup} --no_cstats --insert_sdevs {sdev}".format(min_sup =min_pair_support ,sdev =insert_sdevs ,out =name+'_AA_results/'+name,bed =name+'_AA_CNV_SEEDS.bed',bam = name + '.cs.rmdup.bam')
	print(AA_cmd)
	os.system(AA_cmd)
	######## Updating Amplicon Mapping #################################
	AA_list = os.listdir(name+'_AA_results/')
	for file in AA_list:
		if file.endswith('_graph.txt'):
			if contain_discordant(name+'_AA_results/'+file):
				amplicon_number = file.split('_')[2][8:]
				amplicon_mapping[band].append(amplicon_number)
	##################		Graph Cleaner		##############################
	for amplicon_number in amplicon_mapping[band]:
		clean_cmd = "python2 " + '$PreAA/scripts/graph_cleaner.py' + " -g "+ name + "_AA_results/" + name + "_amplicon"+amplicon_number+"_graph.txt " + "--filter_non_everted --max_hop_size 1000 --max_hop_support 999999"
		print(clean_cmd) 
		os.system(clean_cmd)
		mv_cmd = "mv " + name + "_amplicon"+amplicon_number+"_cleaned_graph.txt " +  "graph_files/."
		print(mv_cmd)
		os.system(mv_cmd)
 	

bulk = parsing_AA_graph(bulk_AA_graph_file)
d = {}
match_count = {}
all_discordant = []
for b in band_list:
	for amplicon_number in amplicon_mapping[b]:
		d[b+'_amplicon'+amplicon_number] = parsing_AA_graph('graph_files/' + cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_graph.txt')
		all_discordant  =  all_discordant + d[b+'_amplicon'+amplicon_number]
		d[b+'_amplicon'+amplicon_number], match_count[b+'_amplicon'+amplicon_number] =  compare_bulk_band(bulk, d[b+'_amplicon'+amplicon_number])

for b in band_list:
	for amplicon_number in amplicon_mapping[b]:
		for i in range(len(d[b+'_amplicon'+amplicon_number])):
			for j in range(len(all_discordant)):
				if compare_discordants(d[b+'_amplicon'+amplicon_number][i], all_discordant[j]):
					d[b+'_amplicon'+amplicon_number][i].count += 1
for b in band_list:
	for amplicon_number in amplicon_mapping[b]:
		with open('graph_files/' + cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_graph.txt' ,'r') as f:
			with open ('filtered_graph_files/' + cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_filtered_graph.txt', 'w') as g:
				for line in f :
					g.write(line)

with open('edge_comparison.txt', 'w') as f:
	for b in band_list:
		for amplicon_number in amplicon_mapping[b]:
			f.write('In band {C}_amplicon{D}, {A} out of {B} are matched to bulk\n'.format(A = match_count[b+'_amplicon'+amplicon_number], B = len(d[b+'_amplicon'+amplicon_number]), C = b,D = amplicon_number))
######################################################################## Finding Path
os.system('mkdir candidate_cycles')
os.system('mkdir candidate_cycles/beds')
os.system('mkdir candidate_cycles/yaml')
os.system('mkdir candidate_cycles/visualization')
os.system('mkdir band_cov')
for b in band_list:
	for amplicon_number in amplicon_mapping[b]:
		find_path_cmd = "python3 {script} -g {graph} --keep_all_LC --remove_short_jumps --runmode isolated --max_length {max_length} --min_length {min_length}".format(script = '$PreAA/scripts/plausible_paths.py', graph = 'filtered_graph_files/' + cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_filtered_graph.txt', max_length=band_size_max[b],min_length = band_size_min[b])
		os.system(find_path_cmd)
		print(find_path_cmd)
		move_cmd = 'mv ' + cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_filtered_candidate_cycles.txt candidate_cycles/.' 
		os.system(move_cmd)
		print(move_cmd)
		generate_cnd_cmd = 'python3 '+ '$PFGE/utils/generate_cnv.py' + ' -i {input} -o {output}'.format(input ='candidate_cycles/'+cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_filtered_candidate_cycles.txt', output = 'candidate_cycles/beds/'+cell_line + '_' + b+'_amplicon'+amplicon_number )
		print(generate_cnd_cmd)
		os.system(generate_cnd_cmd)
quality_report()
######################################################################## visualization
for b in band_list:
	extract_bedgraph_cmd = 'python3 {extract_bedgraph} --bed {bed} --bam {bam} -o {out}'.format(extract_bedgraph = '$CV_SRC/extract_bedgraph.py' , bed = amplicon_bed_file , bam = cell_line + '_' + b+'.cs.rmdup.bam',out = 'band_cov/'+cell_line+'_'+b)
	print(extract_bedgraph_cmd)
	os.system(extract_bedgraph_cmd)
	band_norm_cmd = 'python3 {band_norm} -i {band_cov} -o {out}'.format(band_norm = '$PFGE/utils/band_norm.py' , band_cov ='band_cov/'+cell_line+'_'+b+'_position_coverage.bedgraph', out = 'band_cov/'+cell_line+'_'+b+'_norm.bed' )
	print(band_norm_cmd)
	os.system(band_norm_cmd)

bed_list = os.listdir('candidate_cycles/beds/')
for i in bed_list:
	if i.startswith(cell_line):
		cycle_number = i.split('_')[3].split('.')[0][5:]
		band = i.split('_')[1]
		amplicon_number = i.split('_')[2][8:]
		write_yml('candidate_cycles/yaml/'+cell_line+'_'+band+'_cycle'+cycle_number, 'band_cov/'+cell_line+'_'+band+'_norm.bed')
		cycle_vis_cmd = 'python2 {cycle_viz} --cycles_file {cycle_file} --cycle {cycle} --graph {graph} --ref hg19 --feature_yaml_list {yaml} --label_segs numbers --center_hole 5 --feature_ref_offset 1.5 --noPDF --rotate_to_min -o {output}'.format(
			cycle_viz ='$CV_SRC/CycleViz.py' , cycle_file = 'candidate_cycles/'+cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_filtered_candidate_cycles.txt' , 
			cycle = cycle_number, graph = 'filtered_graph_files/' + cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_filtered_graph.txt',
			yaml = 'candidate_cycles/yaml/'+cell_line+'_'+band+'_cycle'+cycle_number+'.yaml',
			output = 'candidate_cycles/visualization/'+cell_line+'_'+band+'_amplicon'+amplicon_number+'_cycle'+cycle_number)
		print(cycle_vis_cmd) 
		os.system(cycle_vis_cmd)

