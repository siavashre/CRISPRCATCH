# CRISPRCATCH
A tool for analyzing ecDNA structure using data from PFGE amplicon separation. 

## Prerequisites and installation:
CRISPRCATCH has been tested on Ubuntu 20.04.

This tool uses pre-developed tools and the user must install them before using this tool. These tool are:
- [jluebeck/PrepareAA](https://github.com/jluebeck/PrepareAA) is mandatory and bashrc variable `$PreAA` should set to the directory you install PrepareAA. For adding this variable you can change your directory to where you install PrepareAA and run the following commands:
        
        git clone https://github.com/jluebeck/PrepareAA
        cd PrepareAA
        echo export PreAA=$PWD >> ~/.bashrc
        source ~/.bashrc
- The [jluebeck/AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) is mandatory and variable `$AA_SRC` should set in the bashrc. You can find more information about how to install AA and set `$AA_SRC` on the AA Github page.
- The [jluebeck/CycleViz](https://github.com/jluebeck/CycleViz) is mandatory and the variable `$CV_SRC` should set in the bashrc. You can find more information about how install CycleViz and set `$CV_SRC` on CycleViz Github page.
- [samtools](http://www.htslib.org/) (CRISPRCATCH supports versions >= 1.0)
- [bwa mem](https://github.com/lh3/bwa) (BWA MEM should be on the system path after installing so it can be called by PrepareAA)
- [picard](https://github.com/broadinstitute/picard). Variable `$PICARD` should be set in the bashrc file to the directory picard is installed. For doing that, change the directory to where you have picard jar file and run the following commands:

        cd path/to/picard/
        echo export PICARD=$PWD >> ~/.bashrc
        source ~/.bashrc
         
- [intervaltree](https://github.com/chaimleib/intervaltree) (`pip install intervaltree`)

After installing prerequisites, please clone the PFGE repo consider running the following to add a line to your `.bashrc` file:

        cd PFGE/
        echo export PFGE=$PWD >> ~/.bashrc
        source ~/.bashrc
        
## Testing
To ensure that you installed all requiements correctly, you may download these two PFGE CRISPRCATCH test fastq files: [fastq1](https://drive.google.com/file/d/1iYOMtjag3mnZdw5Bqm2cdatE_OwJXFk3/view?usp=sharing) and [fastq2](https://drive.google.com/file/d/1-Vbj6lAsQtQyeXZyT2jZi08HiPT3uUaJ/view?usp=sharing). Then run the following command. If every things installed properly, after ~30 min the pipeline should generate the GBM39 amplicon structure.

`python3 wrapper.py -f1 /path/to/fastq1 -f2 /path/to/fastq2 -b i -sname GBM39 -o /path/to/output/dir -t 1 -r hg19 -bulk test_data/GBM39_AA_graph.txt -lmax 1370 -lmin 1200 -chr chr7 -g_start 55256396 -g_end 55256397`


## Usage
PFGE takes as an input sequence data comes from PFGE experiment for a band. Sequence data in `.fastq` format, expected length of band and outputed best candidates structure for band's structure. 

### Command line arguments
- `-f1` Path to band_r1.fastq file
- `-f2` Path to band_r2.fastq file
- `-b` Analyzed band label name
- `-sname` Analyzed sample name
- `-o` Output folder directory for saving results
- `-t` Number of threads
- `-r` Reference genome version. It can be hg19 or hg38.
- `-bulk` Path to the AA breakpoint graph file for bulk WGS (non-AA)
- `-lmax` Maximum estimated length for structure in terms of kbp. 
- `-lmin` Minimum estimated length for structure in terms of kbp. 
- `-chr` Chromosome of target site. Please write 'chr' + chromosome number.
- `-g_start` CRISPR target site starting position (in bp).
- `-g_end` CRISPR target site ending position (in bp).
- `-bed` path to bed file that spcifying amplicon region. It is not a required argument, PFGE can generate this from AA breakpoint graph file on bulk. 
- `-min_sup` specify the minimum number of reads for calling a breakpoint. Default value is 2. If you have high coverage data (e.g. Novaseq) please set it to 4.
- `-sdv` specify insert_size standard deviation for filtering reads for calling breakpoints. Default value is 8.5. This should be raised if there are an excess of artifactual breakpoints detected by AA.

example command: 

`python3 wrapper.py -f1 SNU16_i_R1.fastq -f2 SNU16_i_R2.fastq -b i -sname SNU16 -o output/ -t 10 -r hg19 -bulk SNU16_AA_amplicon1_graph.txt -lmax 1810 -lmin 1660 -chr chr10 -g_start 123353331 -g_end 123353350 -sdv 8.5 -min_sup 2`

## s_wrapper.py (super wrapper)
If you want to run the pipeline for bunch of bands together (multiple fastq pairs), you can use 's_wrapper.py'.
### Command line arguments
- `-sname` Analyzed Cell-line name
- `-o` Output folder directory for saving results
- `-t` Number of threads
- `-r` Reference genome version. It can be hg19 or hg38.
- `-bulk` path to the AA breakpoint graph file for bulk cell-line
- `-bed`  path to bed file that spcifying amplicon region. It is not required field, PFGE can generate this from AA breakpoint graph file on bulk. 
- `-min_sup` specifying minimum number of reads for calling a breakpoint. Default value is 2. If you have high coverage data please set it to 4.
- `-sdv` specifying insert_size sdv for filtering reads for calling breakpoints. Default value is 8.5.
- `-csv` specifying a csv file containing bands information. Please take a look at 'example.csv' and fill the information as needed. Please use the absolute path for column read1 and read2.
example command: 

`python3 pipeline.py -csv SNU16_i.csv -sname SNU16 -t 10 -r hg19 -o output/ -bulk SNU16_AA_amplicon1_graph.txt -sdv 8.5 -min_sup 2 -bed DNU16.bed`

