import os
from collections import defaultdict
import argparse

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


parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--fastq1", help="fastq file 1", required=True)
parser.add_argument("-f2", "--fastq2", help="fastq file 2", required=True)
parser.add_argument("-b", "--band", help="Band Label", required=True)
parser.add_argument("-sname", "--sname", help="Cell Line name", required=True)
parser.add_argument("-o", "--output", help="Output dir for saving results", required=True)
parser.add_argument("-bed", "--bed", help="bed file describing amplicon region", required=True)
parser.add_argument("-t", "--t", help="Number of thread", required=True)
parser.add_argument("-r", "--ref", help="Reference genome version", required=True)
parser.add_argument("-bulk", "--bulk", help="AA breakpoint graph file for bulk cell line", required=True)
parser.add_argument("-l", "--length", help="Maximum estimated length for this band", required=True)
args = parser.parse_args()
amplicon_mapping = []
output_ans = args.output
amplicon_bed_file = args.bed
generate_cnd_dir = '/nucleus/projects/sraeisid/AA/scripts/generate_cnv.py'
# collapse_cov_dir = '/nucleus/projects/sraeisid/AA/scripts/collapse_cov.pl'
band_norm_dir = '/nucleus/projects/sraeisid/AA/scripts/band_norm.py'

T = 101 #Comparing discordant edges
########################################

band_size = int(args.length)
band = args.band
cell_line = args.sname
os.chdir(output_ans)
name = cell_line + '_' + band
f1 = args.fastq1
f2 = args.fastq2
ref_v = args.ref
Pre_AA_cmd = "python3 {PrepareAA} --ref {ref_v} -t {thread} --use_old_samtools -s {name} --fastqs {f1} {f2} --no_filter --cngain 0 --cnsize_min 0 --cnv_bed {amplicon_bed_file} ".format(
PrepareAA = "$PreAA/PrepareAA.py" , ref_v= ref_v, thread = args.thread, name = name , f1 = f1 , f2 = f2, amplicon_bed_file = amplicon_bed_file)
print(Pre_AA_cmd)
# os.system(Pre_AA_cmd)
os.mkdir(name+'_AA_results')
amplified_amplicon_cmd = '$AA_SRC/amplified_intervals.py --no_cstats --bed {bed} --bam {bam} --ref hg19 --cnsize_min 0 --gain 0 --out {out}'.format(bed=amplicon_bed_file,bam = name + '.cs.rmdup.bam', out= name + '_AA_CNV_SEEDS' )
print(amplified_amplicon_cmd)
# os.system(amplified_amplicon_cmd)
#######		Insert Size Distribution 	########
bam_file = name+'.cs.rmdup.bam'
AA_cmd = "$AA_SRC/AmpliconArchitect.py --out {out} --downsample -1 --bed {bed} --bam {bam} --ref hg19 --pair_support_min 3 --no_cstats --insert_sdevs 8.5".format(out =name+'_AA_results/'+name,bed =name+'_AA_CNV_SEEDS.bed',bam = name + '.cs.rmdup.bam')
print(AA_cmd)
# os.system(AA_cmd)
######## Updating Amplicon Mapping #################################
AA_list = os.listdir(name+'_AA_results/')
for file in AA_list:
	if file.endswith('_graph.txt'):
		if contain_discordant(name+'_AA_results/'+file):
			amplicon_number = file.split('_')[2][8:]
			amplicon_mapping.append(amplicon_number)
##################		Graph Cleaner		##############################
for amplicon_number in amplicon_mapping:
	clean_cmd = "python2 $PreAA/scripts/graph_cleaner.py" + " -g "+ name + "_AA_results/" + name + "_amplicon"+amplicon_number+"_graph.txt " + "--filter_non_everted --max_hop_size 1000 --max_hop_support 999999"
	print(clean_cmd) 
	os.system(clean_cmd)

bulk = parsing_AA_graph(bulk_AA_graph_file)
d = {}
match_count = {}
for amplicon_number in amplicon_mapping:
	d[band+'_amplicon'+amplicon_number] = parsing_AA_graph(cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_graph.txt')
	d[band+'_amplicon'+amplicon_number], match_count[band+'_amplicon'+amplicon_number] =  compare_bulk_band(bulk, d[band+'_amplicon'+amplicon_number])

with open('report.txt', 'w') as f:
	for amplicon_number in amplicon_mapping:
		f.write('In band {C}_amplicon{D}, {A} out of {B} are matched to bulk\n'.format(A = match_count[band+'_amplicon'+amplicon_number], B = len(d[band+'_amplicon'+amplicon_number]), C = band,D = amplicon_number))

######################################################################## Finding Path
os.system('mkdir beds')
for amplicon_number in amplicon_mapping:
	find_path_cmd = "python3 {script} -g {graph} --keep_all_LC --remove_short_jumps --runmode isolated --max_length {max_length}".format(script = '$PreAA/scripts/plausible_paths.py', graph =  cell_line + '_' + b + '_amplicon'+amplicon_number+'_cleaned_graph.txt', max_length=band_size)
	# os.system(find_path_cmd)
	print(find_path_cmd)
	generate_cnd_cmd = 'python3 '+ generate_cnd_dir + ' -i {input} -o {output}'.format(input =cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_candidate_cycles.txt', output = 'beds/'+cell_line + '_' + band+'_amplicon'+amplicon_number )
	print(generate_cnd_cmd)
	# os.system(generate_cnd_cmd)

######################################################################## visualization
extract_bedgraph_cmd = 'python3 {extract_bedgraph} --bed {bed} --bam {bam} -o {out}'.format(extract_bedgraph = '$CV_SRC/extract_bedgraph.py' , bed = amplicon_bed_file , bam = cell_line + '_' + band+'.cs.rmdup.bam',out = cell_line+'_'+band)
print(extract_bedgraph_cmd)
# os.system(extract_bedgraph_cmd)
#Old  band_norm_cmd = 'python3 {band_norm} -i {band_cov} -o {out}'.format(band_norm = band_norm_dir , band_cov ='band_cov/'+cell_line+'_'+b+'.bed', out =  'band_cov/'+cell_line+'_'+b+'_norm.bed')
band_norm_cmd = 'python3 {band_norm} -i {band_cov} -o {out}'.format(band_norm = band_norm_dir , band_cov =cell_line+'_'+band+'_position_coverage.bedgraph', out = cell_line+'_'+band+'_norm.bed' )
print(band_norm_cmd)
# os.system(band_norm_cmd)

bed_list = os.listdir('beds/')
for i in bed_list:
	if i.startswith(cell_line):
		cycle_number = i.split('_')[3].split('.')[0][5:]
		amplicon_number = i.split('_')[2][8:]
		write_yml(cell_line+'_'+band+'_cycle'+cycle_number, cell_line+'_'+band+'_norm.bed')
		cycle_vis_cmd = 'python2 {cycle_viz} --cycles_file {cycle_file} --cycle {cycle} --graph {graph} --ref hg19 --rotate_to_min --feature_yaml_list {yaml} --label_segs numbers --center_hole 5 --feature_ref_offset 1.5 --noPDF -o {output}'.format(
			cycle_viz ='$CV_SRC/CycleViz.py' , cycle_file = cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_candidate_cycles.txt' , 
			cycle = cycle_number, graph =  cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_graph.txt',
			yaml = cell_line+'_'+band+'_cycle'+cycle_number+'.yaml',
			output = cell_line+'_'+band+'_amplicon'+amplicon_number+'_cycle'+cycle_number)
		print(cycle_vis_cmd) 
		# os.system(cycle_vis_cmd)
