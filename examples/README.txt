hmm_cnv.xml is the main configuration file.  It is read in by the executable
analyzer.  The configuration file needs to read a file that contains the input
files and various meta data about that dataset.  In our case it is called
filelist.

Note that filelist contains the sequence of inputfiles that will be analyzed. Be sure that each file match up in the number of SNPs (e.g. same chromosome).

Please unzip the gzipped input files using 

gunzip -v *.gz

before running the example.

Type ./run.sh to run this example of Chr 1 on two breast cancer dilution
series samples as described in our manuscript.
