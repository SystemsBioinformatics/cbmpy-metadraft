# InParanoid-DIAMOND

Accurate and fast ortholog detection.

InParanoid-DIAMOND identifies complex orthologous relationships between protein sequences from different genomes. Through the implementation of DIAMOND in addition to BLAST as the default underlying sequence analysis tool, InParanoid-DIAMOND reduces the runtime of InParanoid by up to 93% while retaining the same high reliability on the orthologs detected. The package is capable of using either DIAMOND or BLAST scores to measure relatedness of proteins, and assigns confidence values for all paralogs in each group. InParanoid is also able to calculate confidence for orthologs using a bootstrap approach.


## Install InParanoid

### Run using Docker container

InParanoid is available as a docker container, including all dependencies needed to run the program. To use the docker container, you need to have docker installed on your local machine. 
See the [Docker website](https://docs.docker.com/engine/install/) for install instructions on different systems. The InParanoid container is managed and maintained on [Dockerhub](https://hub.docker.com/r/sonnhammer/inparanoid) 

To download and test the container (note that you may have to use sudo for docker commands): 
	
	docker pull sonnhammer/inparanoid
	docker run sonnhammer/inparanoid -help
	
The -help command will show you all availabe run options. 

To run InParanoid in the container, you need to mount an input and output directory to it using the -v command. This input-directory has to contain all files you want to analyze, and the program will automatically run all files in 
the directory without the need to use input parameters -f1 -f2 or -input-dir. 
You can append any of the other options listed in the help-function when running the program to the end of the line below. 
Note that text within <> must be replaced with absolute paths to your input and output directories. 

The following line will run InParanoid with default settings on the input files located in the directory specified as <path/to/your/input/files>:
 
	docker run -v <path/to/your/input/files>:/input/ -v <path/where/you/want/output/files>:/output/ sonnhammer/inparanoid
	
#### Run container using Singularity

When root-access is not available on the machine where inparanoid is to be run, for example on a supercomupter cluster, [Singularity](https://sylabs.io/) can be used to run the docker container. 
To run the inparanoid docker container using singularity, use the container in dockerhub tagged "singularity". Your input and output directories have to be mounted to the image in the same way as for docker.
If you have singularity installed, you can fetch and run the container using the following command, where you replace the text within <> with paths to your input/output directories:

	singularity run --bind <path/to/your/input/files>:/input/ --bind <path/where/you/want/output/files>:/output/ docker://sonnhammer/inparanoid:singularity 

### Run using source files

If you do not want to use the Docker container to run InParanoid, you can run the Perl source files directly on your system. This will require you to download the code, and install all necessary dependencies.

1.  Clone this repository to your local machine. 
		
		git clone https://bitbucket.org/sonnhammergroup/inparanoid.git
		
2.  Make sure you have DIAMOND or BLAST installed and available in your system path (see "Dependencies" below on how to install).
3.  To run InParanoid with the default settings, on the example data, run the following command. For the complete set of run options, see the list below.
		
		perl inparanoid.pl -input-dir ./testInput/
		
#### Previous versions

To use an older version of InParanoid, download the code here:

[version 5.0](https://bitbucket.org/sonnhammergroup/inparanoid/get/79d8e39bb243d403ce57b699c1ae104a8a640389.zip )

[version 4.2](https://bitbucket.org/sonnhammergroup/inparanoid4/get/359f8ea484ba.zip)

#### Dependencies

If you want to run InParanoid using the source files, you need to install perl dependencies as well as either DIAMOND or BLAST on your local system. On Ubuntu, follow these commands to install the perl dependencies necessary to run the program:

	apt-get install bioperl libmoose-perl libparallel-forkmanager-perl

##### DIAMOND

In order to run the recommended version of InParanoid, using DIAMOND for sequence alignment, you need to install [DIAMOND](https://github.com/bbuchfink/diamond/wiki/2.-Installation) on your system. In Ubuntu, the binary can be used, and can easily be downloaded using the following commands

	wget http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz
	tar xzf diamond-linux64.tar.gz

Before running InParanoid-DIAMOND, set the path to the diamond binary in your system path, or use the "-diamond-path" option to specify the location of it. Note that DIAMOND will automatically parallelize your sequence alignments, using the threads available on your machine.

##### BLAST

NOTE: To use the recommended default version with DIAMOND, skip this section. 

In case you want to use InParanoid with BLAST, note that legacy [BLAST v2.2.18](http://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.18/blast-2.2.18-x64-linux.tar.gz) must be installed. This can be installed by using the following commands


	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.18/blast-2.2.18-x64-linux.tar.gz 
	tar xvfz blast-2.2.18-x64-linux.tar.gz


## Run InParanoid

When running InParanoid you can specify any of the commandline options listed below. The results of your run will by default appear in the "output" directory. By deafult, this will contain an SQLtable file listing all orthologs detected.

### Required input

InParanoid requires you to input two or more proteome files in fasta format. For further guidance on the format of the input files, please see the example files, EC and SC provided under testInput/.
When running the program on two proteomes, use the -f1 and-f2 options to specify the file names. If you want to run the program on more than two proteomes, use the -input-dir option
to specify a path to a directory containing multiple proteomes in fasta format. This will run InParanoid for all pairs of files in the directory.

### Output

InParnoid can output four different types of files containing the orthologs detected, the SQLtable-file, the stats-file, the html-file and the table-file. 
By default, only the SQLtable file is outputted, but by using the commandline options -out-stats, -out-html and -out-table, these files can be outputted as well.
In addition to the files containing the resulting orthologs, the raw hits from the sequence analysis tool (DIAMOND or BLAST) are stored in the SequenceAnalysisResults-dir under your 
output directory. 

#### SQLtable

The SQLtable file is a tab-separated text file containing the groups of orthologs generated for your search. It lists each protein in a separate row, hence 
having multiple rows for each ortholog group. The file contains the following columns:

	Group-id	Score	Species	Confidence-score	Protein-name
	
#### Stats

The Stats-file is a text-file containing deailed information on the seqences, homologs, matches, orthologs and in-paralogs found in each proteome. The groups of orthologs are listed together with 
its members and scores.

#### HTML

The HTML-file contains similar information as the stats-file but in an html-format, readable by your browser. It contains details on the orthologs and in-paralogs found in each proteome, and lists 
all ortholog groups and its members.

#### Table

The table-file is a tab-separated text file containing all groups of orhologs generated from your seach. This file lists all members of a group in the same row, giving one row per ortholog group
containing one or more proteins from each species. The file contains the following columns

	Group-id	Score	Othologs-from-species-A Confidence-score	Othologs-from-species-B Confidence-score	
	
#### Raw output from sequence analysis tool

The raw output from the sequence analysis tool contains the hits used for the InParanoid analysis. The files are named according to the species run against each other, and four 
files will be generated for each species pair. The files contains the following columns 

	ProteinA	ProteinB	Score	LengthA	LengthB	Length-MatchingRegionA	Length-MatchingRegionB	Total-Length-of-MatchA	Total-Length-of-MatchB		Location-of-Matches
	
#### allPairs

The allPairs file is a tab-separated text file containing all pairs of orthologs generated from the SQLtable files present in the output directory.

	ProteinA	ProteinB


### Command line options

| Commandline option  | Comment |
| ------------- | ------------- |
| - f1  | Fasta file with protein sequences of species A  |
| - f2  | Fasta file with protein sequences of species B  |
| - outgroup  | Fasta file with protein sequences of species C to use as outgroup [Default: no outgroup]  |
| - input-dir  | Directory containing fasta files for multiple species. Will run all vs all. If this option is used, leave -f1 and -f2 empty. Note that InParanoid will run species pairs sequentially,but Diamond will paralellize the sequence search using all available threads. |
| - out-dir  | Specify a directory for the output files [Default: ./output] |
| - seq-tool  | Sequence similarity tool to use. Options: Diamond, Blast [Default: Diamond]  |
| - 2pass  | Run 2-pass approach. Not suitable for Diamond, recommended for Blast [Default: False]  |
| - bootstrap  | Run bootstrapping to estimate confidence of orthologs [Default: False]  |
| - seedscore  | Include calculation of Seed Score to estimate confidence of orthologs [Default: False]  |
| - score-cutoff  | Set bitscore cutoff. Any match below this is ignored [Default: 40] |
| - seq-cutoff  | Set sequence overlap cutoff. Match area should cover at least this much of longer sequence. Match area is the area from start of first segment to end of last segment [Default: 0.5]  |
| - seg-cutoff  | Set segment coverage cutoff. Matching segments must cover this much of thelonger sequence [Default: 0.25]  |
| - outgrp-cutoff  | Set outgroup bitscore cutoff. Outgroup sequencehit must be this many bits stronger to reject best-best hit between A and B [Default: 50]  |
| - conf-cutoff  | Set confidence cutoff. Include in-paralogs with this confidence or better [Default: 0.05]  |
| - grp-cutoff  | Set group overlap cutoff. Merge groups if ortholog in one group has more than this confidence in other group [Default: 0.5]  |
| - grey-zone  | Set grey-zone. This many bits signifies the difference between 2 scores [Default: 0]  |
| - sensitivity  | Set sensitivity mode for Diamond. Options: mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive. [Default: very-sensitive] |
| - matrix  | Specify a matrix to use when running Blast. Options: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, PAM70 [Default: BLOSUM62]  |
| - out-stats  | Output statistics file [Default: False]|
| - out-table  | Output tab-delimited table of orthologs to file [Default: False] |
| - out-sqltable  | Output sqltable file with orthologs [Default: True] |
| - out-html  | Output html file with groups of orthologs [Default: False] |
| -out-allPairs  | Output allPairs file collecting all ortholog pairs from all SQLtable files present in the output directory. [Default: False]  |
| - keep-seqfiles  | Use this option to keep the resulting sequence toolfiles in the working directory. This will let you run InParanoid without re-running the sequence similarity tool. If false, these files will be moved to the output dir when done [Default: False]  |
| - noOverwrite  | Use this option to skip running InParanoid for species-pairs if a resulting SQLtable file is already in the output directory [Default: False]  |
| - diamond-path  | Explicitly state path to Diamond. Can be used if Diamond is in a non-standard location, and not in user PATH [DEFAULT: diamond]  |
|- blast-path  | Explicitly state directory containing blastall and formatdb. Can be used if Blast is in a non-standard location, and not in user PATH.  |
| - cores  | Use to specify the available cores. If DIAMOND is used and this number is higher than twice the -cores-diamond parameter, this number will be split by -cores-diamond to run multiple instances of InParanoid in paralell. If the number is lower, orif only one proteome-pair is run, all cores will be used to run DIAMOND. If BLAST is used, this number will specify the number of paralell InParanoid instances. [Default: using all available cores]  |
| -cores-diamond  | Use to specify the number of cores to use for each DIAMOND run. To optimize performance, please make sure that this number is dividable by the total number of cores used [Default: 4]  |
| - debug  | Activate debug mode [Default: False]  |
| - notimes  | Hide execution times [Default: False] |
| - help/-h  | Show help  |


## References

	Persson E and Sonnhammer E.L.L. (2022)
	InParanoid-DIAMOND: faster orthology analysis with the InParanoid algorithm.
	Bioinformatics. https://doi.org/10.1093/bioinformatics/btac194
	
	Sonnhammer ELL and Ostlund G. (2015)
	InParanoid 8: orthology analysis between 273 proteomes, mostly eukaryotic.
	Nucleic Acids Res. 43:D234-9.


	Berglund AC, Sjolund E, Ostlund G, Sonnhammer ELL (2008)
	InParanoid 6: eukaryotic ortholog clusters with inparalogs"
	Nucleic Acids Res. 36:D263-266


	O'Brien Kevin P, Remm Maido and Sonnhammer Erik L.L (2005). 
	Inparanoid: A Comprehensive Database of Eukaryotic Orthologs
	Nucleic Acids Res. 33:D476-D480


	Remm,M., Storm, C. and E. Sonnhammer (2001). 
	Automatic Clustering of Orthologs and In-paralogs from Pairwise Species Comparisons. 
	J.Mol.Biol. 314:1041-1052


## Contact

If you want to report bugs or if you have any questions, please send an email to <erik.sonnhammer@scilifelab.se> or <emma.persson@scilifelab.se>.


## License

Distributed under the GNU General Public License (GPLv3). See file COPYING


