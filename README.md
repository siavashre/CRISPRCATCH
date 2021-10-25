# CRISPRCATCH
A tool for analyzing ecDNA structure. It works with data comes from PFGE experiment. 

## Prerequisites:
This tool use some pre developed tool and must install them before using this tool. These tool are:
- The [jluebeck/PrepareAA](https://github.com/jluebeck/PrepareAA) is manadtory for PFGE and variable `$PreAA` should set to the directory you install PrepareAA. For adding this variable you can change your directory to where you install PrepareAA and run the following commands:
        
        echo export PreAA=$PWD >> ~/.bashrc
        source ~/.bashrc
- The [jluebeck/AmpliconArchitec](https://github.com/jluebeck/AmpliconArchitect) is manadtory for PFGE and variable `$AA_SRC` should set in the bashrc. You can find more information about how to set `$AA_SRC` on AA Github page.
- The [jluebeck/CycleViz](https://github.com/jluebeck/CycleViz) is manadtory for PFGE and variable `$CV_SRC` should set in the bashrc. You can find more information about how to set `$CV_SRC` on CycleViz Github page.
- [samtools](http://www.htslib.org/) (PFGE supports versions >= 1.0)
- [bwa mem](https://github.com/lh3/bwa)
- [picard](https://github.com/broadinstitute/picard). Variable `$PICARD` should set the directory picard is installed.
- [intervaltree](https://github.com/chaimleib/intervaltree)
## Installation
After installing prerequisites, please clone the PFGE repo consider running the following to add a line to your `.bashrc` file:

        cd PFGE/
        echo export PFGE=$PWD >> ~/.bashrc
        source ~/.bashrc
## Usage
PFGE takes as an input sequence data comes from PFGE experiment for a band. Sequence data in `.fastq` format, expected length of band and outputed best candidates structure for band's structure. 

### Command line arguments
- `-f1` Path to band_r1.fastq file
- `-f2` Path to band_r2.fastq file
- `-b` Analyzed band label name
- `-sname` Analyzed Cell-line name
- `-o` Output folder directory for saving results
- `-t` Number of threads
- `-r` Reference genome version. It can be hg19 or hg38.
- `-bulk` Path to the AA breakpoint graph file for bulk cell-line
- `-lmax` Maximum estimated length for structure in terms of Kbp. 
- `-lmin` Minimum estimated length for structure in terms of Kbp. 
- `-chr` Chromosome of target cite. Please write 'chr' + chromosome number.
- `-g_start` Target cite starting position in bp.
- `-g_end` Target cite ending position in bp.
- `-bed` bed file that spcifying amplicon region. It is not required field, PFGE can generate this from AA breakpoint graph file on bulk. 
- `-min_sup` specifying minimum number of reads for calling a breakpoint. Default value is 2. If you have high coverage data please set it to 4.
- `-sdv` specifying insert_size sdv for filtering reads for calling breakpoints. Default value is 8.5.

example command: 
`python3 wrapper.py -f1 /nucleus/projects/sraeisid/AA/SNU16_run5_high_cov/fastqs/SNU16_i_R1.fastq -f2 /nucleus/projects/sraeisid/AA/SNU16_run5_high_cov/fastqs/SNU16_i_R2.fastq -b i -sname SNU16 -o /nucleus/projects/sraeisid/AA/SNU16_run5_high_cov/test4/ -t 10 -r hg19 -bulk /nucleus/projects/sraeisid/AA/SNU16_run4_clean/ans_Stom1/SNU16_STOMACH_AA_amplicon11_graph.txt -lmax 1810 -lmin 1660 -chr chr10 -g_start 123353331 -g_end 123353350 -sdv 8.5 -min_sup 2`

