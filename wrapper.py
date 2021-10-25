import os
from collections import defaultdict
import argparse
from subprocess import call
import sys

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
						if segments[seg][0] == args.chr:
							if segments[seg][1]-21 < int(args.g_start) < segments[seg][2]+21:
								if segments[seg][1]-21 < int(args.g_end) < segments[seg][2]+21:
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

def quality_report(band):
	header = ['band','cycle_number', 'band_min_length','band_max_length', 'AA_amplicon_id', 'reconstruction_length', 'percent_breakpoints_matched', 'DBI', 'RMSR', 'FILTER','has_cycle', 'includes cut site', 'insert_size mean', 'insert_size stdev']
	with open('report_'+str(band)+'.csv', 'w') as csv_file:
		csvwriter = csv.writer(csv_file)
		csvwriter.writerow(header)
		files = os.listdir(args.output)
		lines = []
		for f in files:
			if f.endswith('candidate_cycles.txt'):
				file_dir = f
				amplicon_number = f.split('_')[2][8:]
				band_max_length = band_size_max
				band_min_length = band_size_min
				lines.extend(parse_cycle(file_dir, band, amplicon_number,band_max_length,band_min_length))
		lines = sorted(lines, key = lambda x:(x[0],x[1]) )
		csvwriter.writerows(lines)


def run_AA(f1,f2,ref_v):
	Pre_AA_cmd = "python3 {PrepareAA} --ref {ref_v} -t {thread} -s {name} --fastqs {f1} {f2} --no_filter --cngain 0 --cnsize_min 0 --cnv_bed {amplicon_bed_file} ".format(
	PrepareAA = "$PreAA/PrepareAA.py" , ref_v= ref_v, thread = args.t, name = name , f1 = f1 , f2 = f2, amplicon_bed_file = amplicon_bed_file)
	print(Pre_AA_cmd)
	call(Pre_AA_cmd,shell=True)
	if not os.path.exists(name+'_AA_results'):
		os.mkdir(name+'_AA_results')
	amplified_amplicon_cmd = '$AA_SRC/amplified_intervals.py --no_cstats --bed {bed} --bam {bam} --ref hg19 --cnsize_min 0 --gain 0 --out {out}'.format(bed=amplicon_bed_file,bam = name + '.cs.rmdup.bam', out= name + '_AA_CNV_SEEDS' )
	print(amplified_amplicon_cmd)
	call(amplified_amplicon_cmd,shell=True)
	bam_file = name+'.cs.rmdup.bam'
	AA_cmd = "$AA_SRC/AmpliconArchitect.py --out {out} --downsample -1 --bed {bed} --bam {bam} --ref hg19 --pair_support_min {min_sup} --no_cstats --insert_sdevs {sdv}".format(sdv = sdv, min_sup = min_sup ,out =name+'_AA_results/'+name,bed =name+'_AA_CNV_SEEDS.bed',bam = name + '.cs.rmdup.bam')
	print(AA_cmd)
	call(AA_cmd,shell=True)
	insert_size_txt = name + '_AA_results/insert_size.txt'
  insert_size_pdf = name + '_AA_results/insert_size.pdf'
  insert_size_cmd = "java -jar $PICARD/picard.jar CollectInsertSizeMetrics I={bam_file} O={insert_size_txt} H={insert_size_pdf} M=0.5".format(bam_file = bam_file , insert_size_txt = insert_size_txt , insert_size_pdf= insert_size_pdf )
  print(insert_size_cmd)
  call(insert_size_cmd, shell=True)
  i_mean, i_std = parse_insert_size(insert_size_txt)
  return i_mean, i_std

def generated_amplicon_bed_file():
	generate_bed_cmd = 'python2 {grah_to_bed} -g {graph} --unmerged'.format(grah_to_bed ='$PreAA/scripts/graph_to_bed.py', graph = args.bulk )
	print(generate_bed_cmd)
	call(generate_bed_cmd,shell=True)

def detecting_amplion_numbers():
	amplicon_mapping = []
	AA_list = os.listdir(name+'_AA_results/')
	for file in AA_list:
		if file.endswith('_graph.txt'):
			if contain_discordant(name+'_AA_results/'+file):
				amplicon_number = file.split('_')[2][8:]
				amplicon_mapping.append(amplicon_number)
	return amplicon_mapping

def run_graph_cleaner():
	for amplicon_number in amplicon_mapping:
		clean_cmd = "python2 $PreAA/scripts/graph_cleaner.py" + " -g "+ name + "_AA_results/" + name + "_amplicon"+amplicon_number+"_graph.txt " + "--filter_non_everted --max_hop_size 1000 --max_hop_support 999999"
		print(clean_cmd) 
		call(clean_cmd)

def compare_bulk_band_report():
	bulk = parsing_AA_graph(args.bulk)
	d = {}
	match_count = {}
	for amplicon_number in amplicon_mapping:
		d[band+'_amplicon'+amplicon_number] = parsing_AA_graph(cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_graph.txt')
		d[band+'_amplicon'+amplicon_number], match_count[band+'_amplicon'+amplicon_number] =  compare_bulk_band(bulk, d[band+'_amplicon'+amplicon_number])
	with open('report.txt', 'w') as f:
		for amplicon_number in amplicon_mapping:
			f.write('In band {C}_amplicon{D}, {A} out of {B} are matched to bulk\n'.format(A = match_count[band+'_amplicon'+amplicon_number], B = len(d[band+'_amplicon'+amplicon_number]), C = band,D = amplicon_number))
	return d , match_count

def run_path_finder():
	if not os.path.exists('beds'):
		os.mkdir('beds')
	for amplicon_number in amplicon_mapping:
		find_path_cmd = "python3 {script} -g {graph} --keep_all_LC --remove_short_jumps --runmode isolated --max_length {max_length} --min_length {min_length}".format(script = '$PreAA/scripts/plausible_paths.py', graph =  cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_graph.txt', max_length=band_size_max, min_length = band_size_min)
		call(find_path_cmd,shell=True)
		print(find_path_cmd)
		generate_cnd_cmd = 'python3 '+ '$PFGE/utils/generate_cnv.py' + ' -i {input} -o {output}'.format(input =cell_line + '_' + band + '_amplicon'+amplicon_number+'_cleaned_candidate_cycles.txt', output = 'beds/'+cell_line + '_' + band+'_amplicon'+amplicon_number )
		print(generate_cnd_cmd)
		call(generate_cnd_cmd,shell=True)

def detect_coverage():
	extract_bedgraph_cmd = 'python3 {extract_bedgraph} --bed {bed} --bam {bam} -o {out}'.format(extract_bedgraph = '$CV_SRC/extract_bedgraph.py' , bed = amplicon_bed_file , bam = cell_line + '_' + band+'.cs.rmdup.bam',out = cell_line+'_'+band)
	print(extract_bedgraph_cmd)
	call(extract_bedgraph_cmd,shell=True)
	band_norm_cmd = 'python3 {band_norm} -i {band_cov} -o {out}'.format(band_norm = '$PFGE/utils/band_norm.py' , band_cov =cell_line+'_'+band+'_position_coverage.bedgraph', out = cell_line+'_'+band+'_norm.bed' )
	print(band_norm_cmd)
	call(band_norm_cmd,shell=True)

def run_visualization():
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
			call(cycle_vis_cmd,shell=True)

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False



parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--fastq1", help="fastq file 1", required=True)
parser.add_argument("-f2", "--fastq2", help="fastq file 2", required=True)
parser.add_argument("-b", "--band", help="Band Label", required=True)
parser.add_argument("-sname", "--sname", help="Cell Line name", required=True)
parser.add_argument("-o", "--output", help="Output dir for saving results", required=True)
parser.add_argument("-t", "--t", help="Number of thread", required=True)
parser.add_argument("-r", "--ref", help="Reference genome version", required=True)
parser.add_argument("-bulk", "--bulk", help="AA breakpoint graph file for bulk cell line", required=True)
parser.add_argument("-lmax", "--lmax", help="Maximum estimated length for this band", required=True)
parser.add_argument("-lmin", "--lmin", help="Minimum estimated length for this band", required=True)
parser.add_argument("-chr", "--chr", help="Chromosome of target cite", required=True)
parser.add_argument("-g_start", "--g_start", help="Target cite starting pos", required=True)
parser.add_argument("-g_end", "--g_end", help="Target cite end pos", required=True)
parser.add_argument("-bed", "--bed", help="bed file describing amplicon region", required=False)
parser.add_argument("-sdv", "--sdv", help="insert_sdev for running AA. Default is 8.5", required=False)
parser.add_argument("-min_sup", "--min_sup", help="Minimum sup pair reads for calling discordant edges. Default is 2", required=False)
args = parser.parse_args()


T = 101 #Comparing discordant edges
band_size_max = int(float(args.lmax))
band_size_min = int(float(args.lmin))
min_sup = 2
sdv = 8.5

if isfloat(args.sdv):
    sdv = float(arg.sdv)
if isfloat(args.min_sup):
    min_sup = float(args.min_sup)
    
band = args.band
cell_line = args.sname
if not os.path.exists(args.output):
	print('Output folder does not exist')
	sys.exit(0) 
os.chdir(args.output)
name = cell_line + '_' + band
if not args.bed:
	generated_amplicon_bed_file()
	amplicon_bed_file = args.bulk[args.bulk.rfind('/')+1:args.bulk.rfind('.')]+'.bed'
else:
	amplicon_bed_file = args.bed
i_mean, i_std = run_AA(args.fastq1,args.fastq2,args.ref)
amplicon_mapping = detecting_amplion_numbers()
run_graph_cleaner()
d, match_count =  compare_bulk_band_report()
run_path_finder()
quality_report(band)
detect_coverage()
run_visualization()

