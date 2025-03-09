# These are the default options for InParanoid 5 using a diamond backend that we will use to generate a user-config widget.

# Internal non-user defined
INPARANOID_SYSTEM_DEFAULTS = {
    "-f1" : (None, True, "Fasta file with protein sequences of species A"),
    "-f2" : (None, True,"Fasta file with protein sequences of species B"),
    "-outgroup" : ("no outgroup", False,"Fasta file with protein sequences of species C to use as outgroup."),
    "-input-dir" : (None, False, "Directory containing fasta files for multiple species. Will run all vs all. If this option is used, leave -f1 and -f2 empty. Note that InParanoid will run species pairs sequentially, but Diamond will paralellize the sequence search using all available threads."), 
    "-out-dir" : ("./output", False, "Specify a directory for the output files."),
    "-seq-tool" : ("Diamond", False, "Sequence similarity tool to use. Options: Diamond, Blast"),
    "-2pass" : ("False", False, "Run 2-pass approach. Not suitable for Diamond, recommended for Blast"),
    "-out-table" : ("True", True, "Output tab-delimited table of orthologs to file, needed for Metadraft."),
    "-out-sqltable" : ("True", False, "Output sqltable file with orthologs."),
    "-noOverwrite" : ("False", False, "Use this option to skip running InParanoid for species-pairs if a resulting SQLtable file is already in the output directory."),
    "-diamond-path" : ("diamond", False, "Explicitly state path to Diamond. Can be used if Diamond is in a non-standard location, and not in user PATH."),
    "-blast-path" : (None, False, "Explicitly state directory containing blastall and formatdb. Can be used if Blast is in a non-standard location, and not in user PATH."),
    "-cores" : ("using all available cores", False, "Use to specify the available cores. If DIAMOND is used and this number is higher than twice the -cores-diamond parameter, this number will be split by -cores-diamond to run multiple instances of InParanoid in paralell. If the number is lower, or if only one proteome-pair is run, all cores will be used to run DIAMOND. If BLAST is used, this number will specify the number of paralell InParanoid instances."),
    "-cores-diamond" : ("4", False, "Use to specify the number of cores to use for each DIAMOND run. To optimize performance, please make sure that this number is dividable by the total number of cores used."),
    "-debug" : ("False", False, "Activate debug mode."),
    "-notimes" : ("False", False, "Hide execution times."),
    "-help" : (None, False, "Show help.")
}

# Exposed user definable options
INPARANOID_USER_DEFAULTS = {
    "-bootstrap" : ("False", False, "Run bootstrapping to estimate confidence of orthologs."),
    "-seedscore" : ("False", False, "Include calculation of Seed Score to estimate confidence of orthologs."),
    "-score-cutoff" : ("40", False, "Set bitscore cutoff. Any match below this is ignored."),
    "-seq-cutoff" : ("0.5", False, "Set sequence overlap cutoff. Match area should cover at least this much of longer sequence. Match area is the area from start of first segment to end of last segment."),
    "-seg-cutoff" : ("0.25", False, "Set segment coverage cutoff. Matching segments must cover this much of the longer sequence."),
    "-outgrp-cutoff" : ("50", False, "Set outgroup bitscore cutoff. Outgroup sequence hit must be this many bits stronger to reject best-best hit between A and B."),
    "-conf-cutoff" : ("0.05", False, "Set confidence cutoff. Include in-paralogs with this confidence or better."),
    "-grp-cutoff" : ("0.5", False, "Set group overlap cutoff. Merge groups if ortholog in one group has more than this confidence in other group."),
    "-grey-zone" : ("0", False, "Set grey-zone. This many bits signifies the difference between 2 scores."),
    "-sensitivity" : ("very-sensitive", False, "Set sensitivity mode for Diamond. Options: mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive."),
    "-matrix" : ("BLOSUM45", False, "Specify a matrix to use when running Blast. Options: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, PAM70."),
    "-out-stats" : ("False", False, "Output statistics file."),
    "-out-html" : ("True", False, "Output html file with groups of orthologs."),
    "-out-allPairs" : ("False", False, "Output allPairs file collecting all ortholog pairs from all SQLtable files present in the output directory."),
    "-keep-seqfiles" : ("False", False, "Use this option to keep the resulting sequence tool files in the working directory. This will let you run InParanoid without re-running the sequence similarity tool. If false, these files will be moved to the output dir when done."),
}
