# PFGE
A tool for analyzing ecDNA structure. It works with data comes from PFGE experiment. 

### Prerequisites:
This tool use some pre developed tool and must install them before using this tool. These tool are:
- The [jluebeck/PrepareAA](https://github.com/jluebeck/PrepareAA) is manadtory for PFGE and variable `$PreAA` should set to the directory you install PrepareAA. For adding this variable you can change your directory to where you install PrepareAA and run the following commands:
        
        echo export PreAA=$PWD >> ~/.bashrc
        source ~/.bashrc
- The [jluebeck/AmpliconArchitec](https://github.com/jluebeck/AmpliconArchitect) is manadtory for PFGE and variable `$AA_SRC` should set in the bashrc. You can find more information about how to set `$AA_SRC` on AA Github page.
- The [jluebeck/CycleViz](https://github.com/jluebeck/CycleViz) is manadtory for PFGE and variable `$CV_SRC` should set in the bashrc. You can find more information about how to set `$CV_SRC` on CycleViz Github page.
- [samtools](http://www.htslib.org/) (PFGE supports versions >= 1.0)
- [bwa mem](https://github.com/lh3/bwa)
