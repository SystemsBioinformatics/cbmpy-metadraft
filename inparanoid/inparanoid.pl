#! /usr/bin/perl
BEGIN { our $start_run = time(); }
###############################################################################
# InParanoid version 5, InParaDiam
# Copyright (C) 2000-2021 Erik Sonnhammer, Kristoffer Forslund, Isabella Pekkari,
# Ann-Charlotte Berglund, Maido Remm, Emma Persson.
#
# This program is provided under the terms of the GNU GENERAL PUBLIC LICENSE Version 3
#
# Make sure that Perl XML libraries are installed!
#
# NOTE: This script requires either DIAMOND (version 2.0.8 or hgher) or blastall (NCBI BLAST)
# version 2.2.16 or higher, that supports compositional score matrix adjustment (-C2 flag)
# Only use legacy BLAST from https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.18/


###############################################################################
# The program calculates orthologs between 2 datasets of proteins
# called A and B. Both datasets should be in multi-fasta file
# - Additionally, it tries to assign paralogous sequences (in-paralogs) to each
#   thus forming paralogous clusters.
# - Confidence of in-paralogs is calculated in relative scale between 0-100%.
#   This confidence value is dependent on how far is given sequence from the
#   seed ortholog of given group
# - Confidence of groups can be calculated with bootstrapping. This is related
#   to score difference between best hit and second best hit.
# - Optionally it can use a species C as outgroup.


###############################################################################
# Main variables for inparanoid:                                              #
###############################################################################

#input to the program
our $fasta_seq_fileA=$ARGV[0]; #Fasta file with protein sequences of species A   #
our $fasta_seq_fileB=$ARGV[1]; #Fasta file with protein sequences of species B   #
our $fasta_seq_fileC; #Fasta file with protein sequences of species C.           #
                      #This file will be used as outgroup in the program         #

# What do you want the program to do?                                            #
my $run_blast = 0;  # Set to 1 to run the BLAST phase                            #
                    # Requires 'blastall', 'formatdb' (NCBI BLAST2)              #
                    # and parser 'blast_parser.pl'                               #
my $two_passes = 0; # Set to 1 to run 2-pass strategy                            #
                    # (strongly recommended for BLAST, NOT for DIAMOND)          #
my $run_inparanoid = 1;                                                          #
my $use_bootstrap = 0;# Use bootstrapping to estimate the confidence of          #
                      # orthologs. Needs additional programs 'seqstat.jar' and   #
                      #'blast2faa.pl'                                            #
my $use_outgroup = 0; # Use proteins from the third genome as an outgroup        #
                      # Reject best-best hit if outgroup sequence is MORE        #
                      # similar to one of the sequences                          #
                      # (by more than $outgroup_cutoff bits)                     #
my $run_diamond=1;    # Set to 1 to run the DIAMOND phase                        #
                      # Requires 'diamond'                                       #
                      # and parser 'diamondParser.pl'                            #
my $sensitivity="very-sensitive"; # Set sensitivity mode for DIAMOND             #
                      # Available options: mid-sensitive, sensitive,             #
                      # more-sensitive, very-sensitive,                          #
                      # ultra-sensitive. Note: very-sensitive is the             #
                      # recommended option                                       #

# Define location of files and programs:                                         #
#$blastall = "blastall -VT"; #Remove -VT for blast version 2.2.12 or earlier     #
my $blastall = "blastall -a1";  #Add -aN to use N processors                     #
my $formatdb = "formatdb";                                                       #
my $seqstat = "seqstat.jar";                                                     #
my $blastParser = "blast_parser.pl";                                             #
my $diamond = "diamond"; # Your path to DIAMOND                                  #
my $diamondParser = "diamondParser.pl";                                          #

my $matrix = "BLOSUM62"; # Reasonable default for comparison of eukaryotes.      #
                         # you can also use BLOSUM45 (for prokaryotes),          #
                         # BLOSUM80 (orthologs within metazoa), PAM70 or PAM30   #

# Output options:                                                                #
my $stats = 0;       # table_stats-format output                                 #
my $table = 0;       # Print tab-delimited table of orthologs to file "table.txt"#
                     # Each orthologous group with all inparalogs is on one line #
my $mysql_table = 1; # Print out sql tables for the web server                   #
                     # Each inparalog is on separate line                        #
my $html = 0;        # HTML-format output                                        #

# Algorithm parameters:
# Default values should work without problems.
# MAKE SURE, however, that the score cutoff here matches what you used for BLAST!#
my $score_cutoff = 40;   # In bits. Any match below this is ignored              #
my $outgroup_cutoff = 50;# In bits. Outgroup sequence hit must be this many bits #
                          # stronger to reject best-best hit between A and B     #
my $conf_cutoff = 0.05;   # Include in-paralogs with this confidence or better   #
my $group_overlap_cutoff = 0.5; # Merge groups if ortholog in one group has more #
                                # than this confidence in other group            #
my $grey_zone = 0;  # This many bits signifies the difference between 2 scores   #
my $show_times = 1; # Show times spent for execution of each part of the program #
                    # (This does not work properly)                              #
my $debug = 0;      # Print debugging messages or not. Levels 0,1,2 and 4 exist  #

my $seq_overlap_cutoff = 0.5; # Match area should cover at least this much of    #
                    # longer sequence. Match area is defined as area from start  #
                    # of first segment to end of last segment, i.e               #
                    # segments 1-10 and 90-100 gives a match length of 100.      #
my $segment_coverage_cutoff = 0.25; 	# Actually matching segments must cover    #
                    # this much of longer sequence. For example, segments 1-10   #
                    # and 90-100 gives a total length of 20.                     #

###############################################################################
# No changes should be required below this line                               #
###############################################################################


#################################################
# Preparing and handling commandline-input
#################################################
use Getopt::Long; #used to capture the command line input in a nice clear way.
use File::Copy; # used to copy input files to workdir + moving output fiels.
use Parallel::ForkManager; # Used to run InParanoid in parallel
no warnings 'experimental::smartmatch';
my $help=0;
my $sequence_tool="Diamond";
my $input_dir;
my $output_dir="./output";
my $keep_seqfiles;
my $blast_path;
my $cores=chomp(my $cores = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
my $cores_diamond=4;
my $allPairsFile=0;
my $helpfile = "help";
my $seedscore = 0;
my $noOverwrite = 0;

# Declaring commandline-input
GetOptions ("f1:s" => \$fasta_seq_fileA,
	    "f2:s" => \$fasta_seq_fileB,
      "outgroup:s" => \$fasta_seq_fileC,
      "input-dir:s" => \$input_dir,
      "out-dir:s" => \$output_dir,
      "seq-tool:s" => \$sequence_tool,
      "2pass!" => \$two_passes,
      "bootstrap!" => \$use_bootstrap,
			"seedscore!" => \$seedscore,
			"no-overwrite!" => \$noOverwrite,
      "matrix:s" => \$matrix,
      "out-stats!" => \$stats,
      "out-table!" => \$table,
      "out-sqltable!" => \$mysql_table,
			"out-html!" => \$html,
			"out-allPairs!" => \$allPairsFile,
			"keep-seqfiles!" => \$keep_seqfiles,
      "score-cutoff:i" => \$score_cutoff,
      "seq-cutoff:f" => \$seq_overlap_cutoff,
      "seg-cutoff:f" => \$segment_coverage_cutoff,
      "outgrp-cutoff:i" => \$outgroup_cutoff,
      "conf-cutoff:f" => \$conf_cutoff,
      "grp-cutoff:f" => \$group_overlap_cutoff,
      "grey-zone:i" => \$grey_zone,
      "sensitivity:s" => \$sensitivity,
      "diamond-path:s" => \$diamond,
      "blast-path:s" => \$blast_path,
			"cores:i" => \$cores,
			"cores-diamond:i" => \$cores_diamond,
      "debug!" => \$debug,
      "times!" => \$show_times,
      "h!" => \$help,
      "help!" => \$help,
)or exit 1;

# Reading content of the help-file
my $document = do {
  local $/ = undef;
  open my $fh, "<", $helpfile
      or die "Could not find the help-file. Has it been accidentally deleted?: $!";
  <$fh>;
};

# Printing help-file if option help, or no input files are given
if ($help){
  print $document;
  exit 1;
}else{
  if ((!$fasta_seq_fileA or !$fasta_seq_fileB) and !$input_dir){
    print "Please use valid commandline options\n",$document;
    exit 1;
  }
}

# Setting use_outgroup to true if outgroup-file is submitted
if($fasta_seq_fileC){
  $use_outgroup=1;
}

# Defining location of blast, if a blast-path is given
if ($blast_path){
  $blastall = "$blast_path/blastall -a1";  #Add -aN to use N processors                     #
  $formatdb = "$blast_path/formatdb";
}

# Checking which sequence similarity tool was selected
if ($sequence_tool eq "Diamond"){
  $run_diamond=1;
  $run_blast=0;
	if (system ("$diamond --version") == -1){
		print "\nERROR: Cant find DIAMOND in system path, please add to your path variable, or specify the location of the DIAMOND binary using -diamond-path \n\n";
		exit 1;
	}
}elsif ($sequence_tool eq "Blast"){
  $run_diamond=0;
  $run_blast=1;
  copy("./matrices/$matrix", "./$matrix");
	if (system ("$blastall >/dev/null 2>&1") == -1){
		print "\nERROR: Cant find BLAST in system path, please add to your path variable, or specify the location of BLAST using -blast-path \n\n";
		exit 1;
	}
	my $formatDBResults=system ("$formatdb >/dev/null 2>&1") ;
	if ($formatDBResults == -1 or $formatDBResults==256){
		print "\nERROR: Cant find BLAST in system path, please add to your path variable, or specify the location of BLAST using -blast-path \n\n";
		exit 1;
	}
	if ($two_passes==0){
		print "\nWARNING: When using Blast for sequence similarity, 2-pass is recommended. activate using parameter -2pass\n\n"
	}
}else{
  print "Please choose one of the following sequence similairty tools: Diamond, Blast\n";
  exit 1;
}

if ($use_bootstrap){
  print "Running with Bootstrap\n" ;
  copy("./matrices/$matrix", "./$matrix");
}

# Verifying diamond sensitivity option selected
if ( !grep { $_ eq $sensitivity} ("mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive") ){
  print "Please coose one of the sensitivity options for Diamond: mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive\n";
  exit 1;
}

# Verifying blast matrix option selected
if ( !grep { $_ eq $matrix} ("BLOSUM62","BLOSUM45", "BLOSUM80", "PAM30", "PAM70") ){
  print "Please coose one of the following matrices: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, PAM70\n";
  exit 1;
}

# Creating output directory if not already exists
if ((-d $output_dir)==False) {
  mkdir $output_dir or die "Error creating directory: $1";
}else{
	print "\nWARNING: writing to existing directory. Note that this may overwrite existing results\n"
}

# Creating directory to put sequence analysis results ( eg. EC-EC files)
if ((-d "$output_dir/SequenceAnalysisResults/")==False) {
  mkdir "$output_dir/SequenceAnalysisResults/" or die "Error creating directory: $1";
}

# If an input directory is given, runs all vs all
my $fasta_seq_fileA_withoutCounter="";
my $fasta_seq_fileB_withoutCounter="";
if ($input_dir){
	if (-e "$input_dir"==False) {
		print("ERROR: Cant find input directory $input_dir\n\n");
		exit 1;
	}
	my @runlist2;
	eval {
	  print "\nWorking with input files in directory: ",$input_dir,"\n";
	  # Generating a list of speciers-pairs to run from all files in the dir
	  my @runlist;
		my @files;
	  opendir(D, $input_dir) || die "Can't open directory: $!\n";
	  while (my $f1 = readdir(D)) {
	    if ($f1 ne "." and $f1 ne ".."){
				if ($f1 ~~ /\s/){
					print "\nERROR: a space was detected in the name of file $f1. Please remove spaces in all file names and try again.\n\n";
					exit 1;
				}
				push @files, $f1;
	    }
	  }
	  closedir(D);

		my @doneFiles;
		opendir(D, $output_dir) || die "Can't open directory: $!\n";
		while (my $f1 = readdir(D)) {
			if ($f1 ne "." and $f1 ne ".." and index($f1,"SQLtable.")>=0 and -s "$output_dir/$f1"){
				push @doneFiles, $f1;
			}
		}
		closedir(D);
		@files2=@files;
		my $filesIndex=1;
		my $len = (scalar @files2)-1;
		foreach $f1 ( @files ) {
			foreach $f2 (@files2[$filesIndex..$len]){
				if ($noOverwrite && (scalar @doneFiles)>0){
					my $valueA="SQLtable.$f1-$f2";
					my $valueB="SQLtable.$f2-$f1";
					if ( grep( /^$valueA$/, @doneFiles ) or grep( /^$valueB$/, @doneFiles ) ) {
  					print "Not running for $f1-$f2, since results are already done\n";
					}else{
						push @runlist, "$f1|$f2";
					}
				}else{
					print "$f1 $f2\n";
					push @runlist, "$f1|$f2";
				}
			}
			$filesIndex+=1;
		}

		# Defining variables for paralell runs
	  @runlist2=@runlist;
	  my $counter=1;
	  my $numberOfRuns= @runlist;
		my $parallelRuns=0;
		if ($run_blast==1){
			$parallelRuns=$cores;
		}elsif ($cores>=($cores_diamond*2)){
			$parallelRuns=int($cores/$cores_diamond);
			$cores=$cores_diamond;
		}else{
			$parallelRuns=0;
		}

		my $pm = Parallel::ForkManager->new($parallelRuns);
		DATA_LOOP:
		# Running in paralell. Note that individual files for each thread needs to be copied and then renamed (for different threads not to overwrite eachothers results)
		for (@runlist){
			my $pid = $pm->start and next DATA_LOOP;
	    print "Running InParanoid for $_ ($counter out of $numberOfRuns) PID=$$","\n";
			my $startTimePair=time;
	    my @files=split('\|', $_);
	    $fasta_seq_fileA="$files[0]:pid$$";
	    $fasta_seq_fileB="$files[1]:pid$$";
			$fasta_seq_fileA_withoutCounter=$files[0];
	    $fasta_seq_fileB_withoutCounter=$files[1];
	    # Always copies inputfiles to current directory
	    copy("$input_dir/$fasta_seq_fileA_withoutCounter","./$fasta_seq_fileA") or die "Copy Failed: $!";
	    copy("$input_dir/$fasta_seq_fileB_withoutCounter","./$fasta_seq_fileB") or die "Copy Failed: $!";
	    runEverything();
	    $counter=$counter+1;

			# Renaming files (removing pids for paralell runs) and removing temp files
			printf ("Running InParanoid for $fasta_seq_fileA - $fasta_seq_fileB took %.2f seconds\n", time - $startTimePair);
			rename("$fasta_seq_fileA-$fasta_seq_fileB","$fasta_seq_fileA_withoutCounter-$fasta_seq_fileB_withoutCounter");
			rename("$fasta_seq_fileB-$fasta_seq_fileA","$fasta_seq_fileB_withoutCounter-$fasta_seq_fileA_withoutCounter" );
			if (-e "$fasta_seq_fileA_withoutCounter-$fasta_seq_fileA_withoutCounter") {
			    unlink "$fasta_seq_fileA-$fasta_seq_fileA";
			} else {
			    rename("$fasta_seq_fileA-$fasta_seq_fileA", "$fasta_seq_fileA_withoutCounter-$fasta_seq_fileA_withoutCounter") || die ( "Error in renaming" );
			}
			if (-e "$fasta_seq_fileB_withoutCounter-$fasta_seq_fileB_withoutCounter") {
			    unlink "$fasta_seq_fileB-$fasta_seq_fileB";
			} else {
			    rename("$fasta_seq_fileB-$fasta_seq_fileB", "$fasta_seq_fileB_withoutCounter-$fasta_seq_fileB_withoutCounter") || die ( "Error in renaming" );
			}
			unlink("./$fasta_seq_fileA");
			unlink("./$fasta_seq_fileB");
			rename("$output_dir/SQLtable.$fasta_seq_fileA-$fasta_seq_fileB", "$output_dir/SQLtable.$fasta_seq_fileA_withoutCounter-$fasta_seq_fileB_withoutCounter") || die ( "Error in renaming" );
			if ($table){
				rename("$output_dir/table.$fasta_seq_fileA-$fasta_seq_fileB", "$output_dir/table.$fasta_seq_fileA_withoutCounter-$fasta_seq_fileB_withoutCounter") || die ( "Error in renaming" );
			}
			if ($stats){
				rename("$output_dir/Stats.$fasta_seq_fileA-$fasta_seq_fileB", "$output_dir/Stats.$fasta_seq_fileA_withoutCounter-$fasta_seq_fileB_withoutCounter") || die ( "Error in renaming" );
			}
			if ($html){
				rename("$output_dir/orthologs.$fasta_seq_fileA-$fasta_seq_fileB.html", "$output_dir/orthologs.$fasta_seq_fileA_withoutCounter-$fasta_seq_fileB_withoutCounter.html") || die ( "Error in renaming" );
			}
			$pm->finish;
	  }
		$pm->wait_all_children;

	  # When everything is done, move seq-analysis files to results dir.
	  if (($keep_seqfiles)==False){
	    for my $r (@runlist2){
	      my @files=split('\|', $r);
	      $fileA=$files[0];
	      $fileB=$files[1];
	      move("./$fileA-$fileA", "$output_dir/SequenceAnalysisResults/$fileA-$fileA");
	      move("./$fileA-$fileB", "$output_dir/SequenceAnalysisResults/$fileA-$fileB");
	      move("./$fileB-$fileB", "$output_dir/SequenceAnalysisResults/$fileB-$fileB");
	      move("./$fileB-$fileA", "$output_dir/SequenceAnalysisResults/$fileB-$fileA");
	    }
	  }
		if ($allPairsFile){
			generateAllPairsFile();
		}
		my ($user,$system,$cuser,$csystem) = times;
		printf ("\nTotal real time: %.2f seconds\n", time - $start_run);
		printf ("Total user time: %.2f seconds\n", $user + $cuser);
		printf ("Total system time: %.2f seconds\n", $system + $csystem);

	} or do {
		print "Something went wrong executing InParanoid. Cleaning up.\n";
	  for my $r (@runlist2){
	    my @files=split('\|', $r);
	    $fileA=$files[0];
	    $fileB=$files[1];
			move("./$fileA-$fileA", "$output_dir/SequenceAnalysisResults/$fileA-$fileA");
			move("./$fileA-$fileB", "$output_dir/SequenceAnalysisResults/$fileA-$fileB");
			move("./$fileB-$fileB", "$output_dir/SequenceAnalysisResults/$fileB-$fileB");
			move("./$fileB-$fileA", "$output_dir/SequenceAnalysisResults/$fileB-$fileA");
			unlink <./$fileA:pid*>;
			unlink <./$fileB:pid*>;
			unlink "$fileA";
			unlink "$fileB";
		}
		exit();
	};
}else{ # If no directory is given, works with one explicit species pair
  # If a path to a file is given, it is copied to the work-dir
	my $copiedFileA=0;
	my $copiedFileB=0;
	my $copiedFileC=0;
	eval {
		# Checking that all files exist AND does not contain space
		if (-e "$fasta_seq_fileA"==False) {
			print("ERROR: Cant find file $fasta_seq_fileA\n");
			exit 1;
		}
		if ((index($fasta_seq_fileA, "/")!=-1 and substr($fasta_seq_fileA, index($fasta_seq_fileA, "/")+1) ~~ /\s/) or $fasta_seq_fileA ~~ /\s/){
			print "\nERROR: a space was detected in the name of file $fasta_seq_fileA. Please remove spaces in all file names and try again.\n\n";
			exit 1;
		}
		if (-e "$fasta_seq_fileB"==False) {
			print("ERROR: Cant find file $fasta_seq_fileB\n");
			exit 1;
		}
		if ((index($fasta_seq_fileB, "/")!=-1 and substr($fasta_seq_fileB, index($fasta_seq_fileB, "/")+1) ~~ /\s/) or $fasta_seq_fileB ~~ /\s/){
			print "\nERROR: a space was detected in the name of file $fasta_seq_fileB. Please remove spaces in all file names and try again.\n\n";
			exit 1;
		}
		if ($fasta_seq_fileC and -e "$fasta_seq_fileC"==False) {
			print("ERROR: Cant find file $fasta_seq_fileC\n");
			exit 1;
		}
		if ($fasta_seq_fileC and ((index($fasta_seq_fileC, "/")!=-1 and substr($fasta_seq_fileC, index($fasta_seq_fileC, "/")+1) ~~ /\s/) or $fasta_seq_fileC ~~ /\s/)){
			print "\nERROR: a space was detected in the name of file $fasta_seq_fileC. Please remove spaces in all file names and try again.\n\n";
			exit 1;
		}
		# Moving files to wokdir if necessary
	  if (index($fasta_seq_fileA, "/") != -1){
	    my @filename=split('\/',$fasta_seq_fileA);
	    copy("$fasta_seq_fileA","./@filename[-1]") or die print "\nERROR: Failed to copy $fasta_seq_fileA to working directory \n\n";
	    $copiedFileA=1;
	    $fasta_seq_fileA=@filename[-1];
	  }
		$fasta_seq_fileA_withoutCounter=$fasta_seq_fileA;
	  if (index($fasta_seq_fileB, "/") != -1){
	    my @filename=split('\/',$fasta_seq_fileB);
	    copy("$fasta_seq_fileB","./@filename[-1]") or die print "\nERROR: Failed to copy $fasta_seq_fileB to working directory \n\n";
	    $copiedFileB=1;
	    $fasta_seq_fileB=@filename[-1];
	  }
		$fasta_seq_fileB_withoutCounter=$fasta_seq_fileB;
	  if ($fasta_seq_fileC and index($fasta_seq_fileC, "/") != -1){
	    my @filename=split('\/',$fasta_seq_fileC);
	    copy("$fasta_seq_fileB","./@filename[-1]") or die print "\nERROR: Please make sure that you are not using the same file for outgroup as for one of the input files \n\n";
	    $copiedFileC=1;
	    $fasta_seq_fileC=@filename[-1];
	  }

	  runEverything();

	  if ($copiedFileA){
	    unlink("./$fasta_seq_fileA");
	  }
	  if ($copiedFileB){
	    unlink("./$fasta_seq_fileB");
	  }
	  if ($copiedFileC){
	    unlink("./$fasta_seq_fileC");
	  }
	  if (($keep_seqfiles)==False){
	    move("./$fasta_seq_fileA-$fasta_seq_fileA", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileA-$fasta_seq_fileA");
	    move("./$fasta_seq_fileA-$fasta_seq_fileB", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileA-$fasta_seq_fileB");
	    move("./$fasta_seq_fileB-$fasta_seq_fileB", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileB-$fasta_seq_fileB");
	    move("./$fasta_seq_fileB-$fasta_seq_fileA", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileB-$fasta_seq_fileA");
	    if ($fasta_seq_fileC){
	      move("./$fasta_seq_fileA-$fasta_seq_fileC", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileA-$fasta_seq_fileC");
	      move("./$fasta_seq_fileB-$fasta_seq_fileC", "$output_dir/SequenceAnalysisResults/$fasta_seq_fileB-$fasta_seq_fileC");
	    }
	  }
		if ($allPairsFile){
			generateAllPairsFile();
		}
		my ($user,$system,$cuser,$csystem) = times;
		printf ("\nTotal real time: %.2f seconds\n", time - $start_run);
		printf ("Total user time: %.2f seconds\n", $user + $cuser);
		printf ("Total system time: %.2f seconds\n", $system + $csystem);
	} or do {
		print "Something went wrong executing InParanoid. Cleaning up.\n";
		if ($copiedFileA){
	    unlink("./$fasta_seq_fileA");
	  }
	  if ($copiedFileB){
	    unlink("./$fasta_seq_fileB");
	  }
		unlink "$fasta_seq_fileA-$fasta_seq_fileA";
		unlink "$fasta_seq_fileB-$fasta_seq_fileB";
		unlink "$fasta_seq_fileA-$fasta_seq_fileB";
		unlink "$fasta_seq_fileB-$fasta_seq_fileA";
		if ($fasta_seq_fileC){
			if ($copiedFileC){
				unlink("./$fasta_seq_fileC");
			}
			unlink $fasta_seq_fileA-$fasta_seq_fileC;
			unlink $fasta_seq_fileB-$fasta_seq_fileC;
		}

	};
}


sub generateAllPairsFile{ #make an all-pairs file from every SQLtable file in the output dir
	my $totalPairs=0;
	open(FH, '>', "allPairs") or die $!;
	opendir(D, $output_dir) || die "Can't open directory: $!\n";
	while (my $f1 = readdir(D)) {
		if (index($f1, "SQLtable") != -1) {
			my $protSpecies;
			my $groups;
			open my $info, "$output_dir/$f1" or die "Could not open $f1: $!";
			while( my $line = <$info>)  {
				my @tokens = split /\t/, $line;
				my $prot=$tokens[4];
				$prot=~ s/^\s+|\s+$//g;
				$protSpecies{$prot}=$tokens[2];
				if (exists $groups{$tokens[0]}) {
					push @{ $groups{$tokens[0]} }, $prot;
				} else {
					$groups{$tokens[0]}=[$prot]
				}
			}
			close $info;
			foreach $key (keys %groups)
			{
			  my @value = @{ $groups{$key} };
				my @added=[];
				foreach $protA (@value){
					foreach $protB (@value){
						if (($protSpecies{$protA} ne $protSpecies{$protB}) && (! ("$protA|$protB" ~~ @added))){
							push (@added, "$protB|$protA");
							print FH "$protA\t$protB\n";
							$totalPairs+=1;
						}

					}
				}
			}
			undef %protSpecies;
			undef %groups;
		}
	}
	closedir(D);
	close(FH);
	print "\nTotal pairs: $totalPairs\n";
	move("./allPairs", "$output_dir/allPairs");
}

#################################################
# Declaring variables to use later
#################################################
sub runEverything { # Putting EVERYTHING in the same subroutine (including
                    # other subrouines, since everyone is using the same global
                    # variables and nothing works otherwise) to be able to run
                    # it for all vs all fasta-files in a directory.
$ENV{CLASSPATH} = "./$seqstat" if ($use_bootstrap);
if ((!$run_blast) and (!$run_inparanoid) and (!$run_diamond)){
    print STDERR "run_blast or run_inparanoid has to be set!\n";
    exit 1;
}

my $seqSimilarity_outputAB = $fasta_seq_fileA . "-" . $fasta_seq_fileB;
my $seqSimilarity_outputBA = $fasta_seq_fileB . "-" . $fasta_seq_fileA;
my $seqSimilarity_outputAA = $fasta_seq_fileA . "-" . $fasta_seq_fileA;
my $seqSimilarity_outputBB = $fasta_seq_fileB . "-" . $fasta_seq_fileB;
my $seqSimilarity_outputAB_alternative = $fasta_seq_fileA_withoutCounter . "-" . $fasta_seq_fileB_withoutCounter;
my $seqSimilarity_outputBA_alternative = $fasta_seq_fileB_withoutCounter . "-" . $fasta_seq_fileA_withoutCounter;
my $seqSimilarity_outputAA_alternative = $fasta_seq_fileA_withoutCounter . "-" . $fasta_seq_fileA_withoutCounter;
my $seqSimilarity_outputBB_alternative = $fasta_seq_fileB_withoutCounter . "-" . $fasta_seq_fileB_withoutCounter;

if ($use_outgroup){
    $seqSimilarity_outputAC = $fasta_seq_fileA . "-" . $fasta_seq_fileC;
    $seqSimilarity_outputBC = $fasta_seq_fileB . "-" . $fasta_seq_fileC;
}
my %idA;        # Name -> ID combinations for species 1
my %idB;        # Name -> ID combinations for species 2
my @nameA;      # ID -> Name combinations for species 1
my @nameB;      # ID -> Name combinations for species 2
my @nameC;
my %scoreAB;    # Hashes with pairwise BLAST scores (in bits)
my %scoreBA;
my %scoreAA;
my %scoreBB;
my @hitnAB;     # 1-D arrays that keep the number of pairwise hits
my @hitnBA;
my @hitnAA;
my @hitnBB;
my @hitAB;      # 2-D arrays that keep the actual matching IDs
my @hitBA;
my @hitAA;
my @hitBB;
my @besthitAB;  # IDs of best hits in other species (may contain more than one ID)
my @besthitBA;  # IDs of best hits in other species (may contain more than one ID)
my @bestscoreAB; # best match A -> B
my @bestscoreBA; # best match B -> A
my @ortoA;      # IDs of ortholog candidates from species A
my @ortoB;      # IDs of ortholog candidates from species B
my @ortoS;      # Scores between ortoA and ortoB pairs
my @paralogsA;  # List of paralog IDs in given cluster
my @paralogsB;  # List of paralog IDs in given cluster
my @confPA;     # Confidence values for A paralogs
my @confPB;     # Confidence values for B paralogs
my @confA;      # Confidence values for orthologous groups
my @confB;      # Confidence values for orthologous groups
my $prev_time = 0;
my $o;


if ($stats && $run_inparanoid){
    $statsfile = "Stats." . $fasta_seq_fileA . "-" . $fasta_seq_fileB;
    open OUTPUT, ">$statsfile" or warn "Could not write to Stats output file $filename\n";
}

#################################################
# Assign ID numbers for species A
#################################################
open A, "$fasta_seq_fileA" or die "File A with sequences in FASTA format is missing
Usage $0 <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <FASTAFILE with sequences of species C>\n";
$id = 0;
while (<A>){
    if(/^\>/){
	++$id;
	chomp;
	s/\>//;
	@tmp = split(/\s+/);
	$name = $tmp[0];
	$idA{$name} = int($id);
	$nameA[$id] = $name;
    }
}
close A;
$A = $id;
print "$A sequences in file $fasta_seq_fileA\n";

if ($stats && $run_inparanoid){
    print OUTPUT "$A sequences in file $fasta_seq_fileA_withoutCounter\n";
}

#{
#################################################
# Assign ID numbers for species B
#################################################
open B, "$fasta_seq_fileB" or die "File B with sequences in FASTA format is missing
Usage $0 <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <FASTAFILE with sequences of species C>\n";
$id = 0;
while (<B>){
	if(/^\>/){
    ++$id;
    chomp;
    s/\>//;
    @tmp = split(/\s+/);
    $name = $tmp[0];
    $idB{$name} = int($id);
    $nameB[$id] = $name;
	}
}
$B = $id;
print "$B sequences in file $fasta_seq_fileB\n";
close B;

if ($stats && $run_inparanoid){
	print OUTPUT "$B sequences in file $fasta_seq_fileB_withoutCounter\n";
}
#}

#################################################
# Assign ID numbers for species C (outgroup)
#################################################
if ($use_outgroup){
	open C, "$fasta_seq_fileC" or die "File C with sequences in FASTA format is missing
  Usage $0 <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <FASTAFILE with sequences of species C>\n";
  $id = 0;
  while (<C>){
		if(/^\>/){
	    ++$id;
	    chomp;
	    s/\>//;
	    @tmp = split(/\s+/);
	    $name = $tmp[0];
	    $idC{$name} = int($id);
	    $nameC[$id] = $name;
		}
  }
  $C = $id;
  print "$C sequences in file $fasta_seq_fileC\n";
  close C;
  if ($stats && $run_inparanoid){
		print OUTPUT "$C sequences in file $fasta_seq_fileC\n";
  }
}

if ($show_times){
    ($user_time,,,) = times;
    printf ("Indexing sequences took %.2f seconds\n", ($user_time - $prev_time));
    $prev_time = $user_time;
}

#################################################
# Run BLAST if not done already
#################################################
if ($run_blast) {

    # Run blast only if the files do not already exist.
    print "Starting to run BLAST ... \n";
    do_blast ($fasta_seq_fileA, $fasta_seq_fileA, $A, $A, $seqSimilarity_outputAA) if !-e $seqSimilarity_outputAA && !-e $seqSimilarity_outputAA_alternative;
    do_blast ($fasta_seq_fileB, $fasta_seq_fileB, $B, $B, $seqSimilarity_outputBB) if !-e $seqSimilarity_outputBB && !-e $seqSimilarity_outputBB_alternative;
    do_blast ($fasta_seq_fileA, $fasta_seq_fileB, $B, $B, $seqSimilarity_outputAB) if !-e $seqSimilarity_outputAB && !-e $seqSimilarity_outputAB_alternative;
    do_blast ($fasta_seq_fileB, $fasta_seq_fileA, $A, $A, $seqSimilarity_outputBA) if !-e $seqSimilarity_outputBA && !-e $seqSimilarity_outputBA_alternative;

    if ($use_outgroup){
			do_blast ($fasta_seq_fileA, $fasta_seq_fileC, $A, $C, $seqSimilarity_outputAC) if !-e $seqSimilarity_outputAC;
			do_blast ($fasta_seq_fileB, $fasta_seq_fileC, $B, $C, $seqSimilarity_outputBC) if !-e $seqSimilarity_outputBC;
    }

    if ($show_times){
			($user_time,,,) = times;
			printf ("BLAST searches took %.2f seconds\n", ($user_time - $prev_time));
			$prev_time = $user_time;
    }
    print STDERR "Done BLAST searches. ";

}

#####################################################################
########## Run Diamond ##############################################
#####################################################################

if ($run_diamond) {
	do_diamond ($fasta_seq_fileA, $fasta_seq_fileA, $A, $A, $seqSimilarity_outputAA) if !-e $seqSimilarity_outputAA && !-e $seqSimilarity_outputAA_alternative;
	do_diamond ($fasta_seq_fileB, $fasta_seq_fileB, $B, $B, $seqSimilarity_outputBB) if !-e $seqSimilarity_outputBB && !-e $seqSimilarity_outputBB_alternative;
	do_diamond ($fasta_seq_fileA, $fasta_seq_fileB, $B, $B, $seqSimilarity_outputAB) if !-e $seqSimilarity_outputAB && !-e $seqSimilarity_outputAB_alternative;
	do_diamond ($fasta_seq_fileB, $fasta_seq_fileA, $A, $A, $seqSimilarity_outputBA) if !-e $seqSimilarity_outputBA && !-e $seqSimilarity_outputBA_alternative;

  if ($use_outgroup){
		do_diamond ($fasta_seq_fileA, $fasta_seq_fileC, $A, $C, $seqSimilarity_outputAC) if !-e $seqSimilarity_outputAC;
  	do_diamond ($fasta_seq_fileB, $fasta_seq_fileC, $B, $C, $seqSimilarity_outputBC) if !-e $seqSimilarity_outputBC;
  }
  if ($show_times){
		($user_time,,,) = times;
		printf ("DIAMOND searches took %.2f seconds\n", ($user_time - $prev_time));
		$prev_time = $user_time;
  }
  print STDERR "Done with DIAMOND searches. ";
}

if ($run_inparanoid){
    print STDERR "Starting ortholog detection...\n";
#################################################
# Read in best hits from blast output file AB
#################################################
    $count = 0;
		if (-e "$seqSimilarity_outputAB"){
			open AB, "$seqSimilarity_outputAB" or die "Blast output file A->B is missing\n";
		}elsif(-e "$seqSimilarity_outputAB_alternative"){
			open AB, "$seqSimilarity_outputAB_alternative" or die "Blast output file A->B is missing\n";
		}else{
			die "Blast output file A->B is missing\n";
		}
    $old_idQ = 0;
    while (<AB>){
			chomp;
			@Fld = split(/\s+/);    # Get query, match and score
			if( scalar @Fld < 9 ){
	    	if($Fld[0]=~/done/){
					print STDERR "AB ok\n";
	    	}
	    next;
			}
			$q = $Fld[0];
			$m = $Fld[1];
			$idQ = $idA{$q}; # ID of query sequence
			$idM = $idB{$m}; # ID of match sequence
			$score = $Fld[2];
			next if (!overlap_test(@Fld));
			# Score must be equal to or above cut-off
			next if ($score < $score_cutoff);

			if(!$count || $q ne $oldq){
		    print "Match $m, score $score, ID for $q is missing\n" if ($debug == 2 and !(exists($idA{$q})));
		    $hitnAB[$idA{$oldq}] = $hit if($count); # Record number of hits for previous query
		    $hit = 0;
		    ++$count;
		    $oldq = $q;
			}
			++$hit;
			$hitAB[$idQ][$hit] = int($idM);
			$scoreAB{"$idQ:$idM"} = $score;
			$scoreBA{"$idM:$idQ"} = $score_cutoff; # Initialize mutual hit score - sometimes this is below score_cutoff
			$old_idQ = $idQ;
    }
    $hitnAB[$idQ] = $hit; # For the last query
    close AB;
    if ($stats){
			print OUTPUT "$count sequences $fasta_seq_fileA_withoutCounter have homologs in dataset $fasta_seq_fileB_withoutCounter\n";
    }
#################################################
# Read in best hits from blast output file BA
#################################################
    $count = 0;
		if (-e "$seqSimilarity_outputBA"){
			open BA, "$seqSimilarity_outputBA" or die "Blast output file B->A is missing\n";
		}elsif(-e "$seqSimilarity_outputBA_alternative"){
			open BA, "$seqSimilarity_outputBA_alternative" or die "Blast output file B->A is missing\n";
		}else{
			die "Blast output file B->A is missing\n";
		}
    $old_idQ = 0;
    while (<BA>){
			chomp;
			@Fld = split(/\s+/);    # Get query, match and score
			if( scalar @Fld < 9 ){
	    	if($Fld[0]=~/done/){
					print STDERR "BA ok\n";
	    	}
	    	next;
			}
			$q = $Fld[0];
			$m = $Fld[1];
			$idQ = $idB{$q};
			$idM = $idA{$m};
			$score = $Fld[2];

			next if (!overlap_test(@Fld));

			next if ($score < $score_cutoff);

			if(!$count || $q ne $oldq){
		    print "ID for $q is missing\n" if ($debug == 2 and (!exists($idB{$q})));
		    $hitnBA[$idB{$oldq}] = $hit if($count); # Record number of hits for previous query
		    $hit = 0;
		    ++$count;
		    $oldq = $q;
			}
			++$hit;
			$hitBA[$idQ][$hit] = int($idM);
			$scoreBA{"$idQ:$idM"} = $score;
			$scoreAB{"$idM:$idQ"} = $score_cutoff if ($scoreAB{"$idM:$idQ"} < $score_cutoff); # Initialize missing scores
			$old_idQ = $idQ;
    }
    $hitnBA[$idQ] = $hit; # For the last query
    close BA;
    if ($stats){
			print OUTPUT "$count sequences $fasta_seq_fileB_withoutCounter have homologs in dataset $fasta_seq_fileA_withoutCounter\n";
    }
##################### Equalize AB scores and BA scores ##########################

###################################################################################################################################### Modification by Isabella 1

	# I removed the time consuming all vs all search and equalize scores for all pairs where there was a hit
    foreach my $key (keys %scoreAB) {
			my ($a, $b) = split(':', $key);
			my $key2 = $b . ':' . $a;
			# If debugg mod is 5 and the scores A-B and B-A are unequal
	   	# the names of the two sequences and their scores are printed
	    if ($scoreAB{$key} != $scoreBA{$key2}){
				printf ("%-20s\t%-20s\t%d\t%d\n",$nameA[$a], $nameB[$b], $scoreAB{$key}, $scoreBA{$key2}) if ($debug == 5);
	    }
			# Set score AB and score BA to the mean of scores AB and BA.
	  	# The final score is saved as an integer so .5 needs to be added to avoid rounding errors
	    $scoreAB{$key} = $scoreBA{$key2} = int(($scoreAB{$key} + $scoreBA{$key2})/2.0 +.5);
		}
####################################################################################################################################### End modification by Isabella 1

##################### Re-sort hits, besthits and bestscore #######################
    for $idA(1..$A){
			next if (!($hitnAB[$idA]));
			for $hit (1..($hitnAB[$idA]-1)){ # Sort hits by score
	    	while($scoreAB{"$idA:$hitAB[$idA][$hit]"} < $scoreAB{"$idA:$hitAB[$idA][$hit+1]"}){
					$tmp = $hitAB[$idA][$hit];
					$hitAB[$idA][$hit] = $hitAB[$idA][$hit+1];
					$hitAB[$idA][$hit+1] = $tmp;
					--$hit if ($hit > 1);
	    	}
			}
			$bestscore = $bestscoreAB[$idA] = $scoreAB{"$idA:$hitAB[$idA][1]"};
			$besthitAB[$idA] = $hitAB[$idA][1];
			for $hit (2..$hitnAB[$idA]){
	    	if ($bestscore - $scoreAB{"$idA:$hitAB[$idA][$hit]"} <= $grey_zone){
					$besthitAB[$idA] .= " $hitAB[$idA][$hit]";
	    	}
	    	else {
					last;
	    	}
			}
			undef $is_besthitAB[$idA]; # Create index that we can check later
			grep (vec($is_besthitAB[$idA],$_,1) = 1, split(/ /,$besthitAB[$idA]));
    }

    for $idB(1..$B){
			next if (!($hitnBA[$idB]));
			for $hit (1..($hitnBA[$idB]-1)){
				# Sort hits by score
	    	while($scoreBA{"$idB:$hitBA[$idB][$hit]"} < $scoreBA{"$idB:$hitBA[$idB][$hit+1]"}){
					$tmp = $hitBA[$idB][$hit];
					$hitBA[$idB][$hit] = $hitBA[$idB][$hit+1];
					$hitBA[$idB][$hit+1] = $tmp;
					--$hit if ($hit > 1);
	    	}
			}
			$bestscore = $bestscoreBA[$idB] = $scoreBA{"$idB:$hitBA[$idB][1]"};
			$besthitBA[$idB] = $hitBA[$idB][1];
			for $hit (2..$hitnBA[$idB]){
	    	if ($bestscore - $scoreBA{"$idB:$hitBA[$idB][$hit]"} <= $grey_zone){
					$besthitBA[$idB] .= " $hitBA[$idB][$hit]";
	    	} else {last;}
			}
			undef $is_besthitBA[$idB]; # Create index that we can check later
			grep (vec($is_besthitBA[$idB],$_,1) = 1, split(/ /,$besthitBA[$idB]));
    }

    if ($show_times){
			($user_time,,,) = times;
			printf ("Reading and sorting homologs took %.2f seconds\n", ($user_time - $prev_time));
			$prev_time = $user_time;
    }

######################################################
# Now find orthologs:
######################################################
    $o = 0;
    for $i(1..$A){   # For each ID in file A
			if (defined $besthitAB[$i]){
	    	@besthits = split(/ /,$besthitAB[$i]);
	    	for $hit (@besthits){
					if (vec($is_besthitBA[$hit],$i,1)){
				    ++$o;
				    $ortoA[$o] = $i;
				    $ortoB[$o] = $hit;
				    $ortoS[$o] = $scoreAB{"$i:$hit"}; # Should be equal both ways
			    	print "Accept! " if ($debug == 2);
					} else {print "        " if ($debug == 2);}
					printf ("%-20s\t%d\t%-20s\t", $nameA[$i], $bestscoreAB[$i], $nameB[$hit]) if ($debug == 2);
					print "$bestscoreBA[$hit]\t$besthitBA[$hit]\n" if ($debug == 2);
	    	}
			}
    }
    print "$o ortholog candidates detected\n" if ($debug);
#####################################################
# Sort orthologs by ID and then by score:
#####################################################

####################################################################################################### Modification by Isabella 2

    # Removed time consuiming bubble sort. Created an index array and sort that according to id and score.
    # The put all clusters on the right place.
    # Create an array used to store the position each element shall have in the final array
    # The elements are initialized with the position numbers
    my @position_index_array = (1..$o);
    # Sort the position list according to id
    my @id_sorted_position_list = sort { ($ortoA[$a]+$ortoB[$a]) <=> ($ortoA[$b] + $ortoB[$b]) } @position_index_array;
    # Sort the list according to score
    my @score_id_sorted_position_list = sort { $ortoS[$b] <=> $ortoS[$a] } @id_sorted_position_list;
    # Create new arrays for the sorted information
    my @new_ortoA;
    my @new_ortoB;
    my @new_orthoS;
   # Add the information to the new arrays in the orer specifeid by the index array
	 for (my $index_in_list = 0; $index_in_list < scalar @score_id_sorted_position_list; $index_in_list++) {
		 	my $old_index = $score_id_sorted_position_list[$index_in_list];
			$new_ortoA[$index_in_list + 1] = $ortoA[$old_index];
			$new_ortoB[$index_in_list + 1] = $ortoB[$old_index];
			$new_ortoS[$index_in_list + 1] = $ortoS[$old_index];
	 }
	 @ortoA = @new_ortoA;
	 @ortoB = @new_ortoB;
	 @ortoS = @new_ortoS;

###################################################################################################### End modification bt Isabella 2

    @all_ortologsA = ();
    @all_ortologsB = ();
    for $i(1..$o){
			push(@all_ortologsA,$ortoA[$i]); # List of all orthologs
			push(@all_ortologsB,$ortoB[$i]);
    }
    undef $is_ortologA; # Create index that we can check later
    undef $is_ortologB;
    grep (vec($is_ortologA,$_,1) = 1, @all_ortologsA);
    grep (vec($is_ortologB,$_,1) = 1, @all_ortologsB);

    if ($show_times){
			($user_time,,,) = times;
			printf ("Finding and sorting orthologs took %.2f seconds\n", ($user_time - $prev_time));
			$prev_time = $user_time;
    }
#################################################
# Read in best hits from blast output file AC
#################################################
    if ($use_outgroup){
			$count = 0;
			open AC, "$seqSimilarity_outputAC" or die "Blast output file A->C is missing\n";
			while (<AC>){
	    	chomp;
	    	@Fld = split(/\s+/);    # Get query, match and score
        if( scalar @Fld < 9 ){
	      	if($Fld[0]=~/done/){
		   			print STDERR "AC ok\n";
   	      }
	      next;
        }
		    $q = $Fld[0];
		    $m = $Fld[1];
		    $idQ = $idA{$q};
		    $idM = $idC{$m};
		    $score = $Fld[2];
		    next unless (vec($is_ortologA,$idQ,1));
		    next if (!overlap_test(@Fld));
		    next if ($score < $score_cutoff);
		    next if ($count and ($q eq $oldq));
		    # Only comes here if this is the best hit:
		    $besthitAC[$idQ] = $idM;
		    $bestscoreAC[$idQ] = $score;
		    $oldq = $q;
		    ++$count;
			}
			close AC;
#################################################
# Read in best hits from blast output file BC
#################################################
			$count = 0;
			open BC, "$seqSimilarity_outputBC" or die "Blast output file B->C is missing\n";
			while (<BC>){
	    	chomp;
	    	@Fld = split(/\s+/);    # Get query, match and score
        if( scalar @Fld < 9 ){
	       if($Fld[0]=~/done/){
		   		print STDERR "BC ok\n";
   	     }
	       next;
        }
		    $q = $Fld[0];
		    $m = $Fld[1];
		    $idQ = $idB{$q};
		    $idM = $idC{$m};
		    $score = $Fld[2];
		    next unless (vec($is_ortologB,$idQ,1));
		    next if (!overlap_test(@Fld));
		    next if ($score < $score_cutoff);
		    next if ($count and ($q eq $oldq));
		    # Only comes here if this is the best hit:
		    $besthitBC[$idQ] = $idM;
		    $bestscoreBC[$idQ] = $score;
		    $oldq = $q;
		    ++$count;
			}
			close BC;
################################
# Detect rooting problems
################################
		$rejected = 0;
		@del = ();
		$file = "rejected_sequences." . $fasta_seq_fileC;
		open OUTGR, ">$file";
		for $i (1..$o){
	    $diff1 = $diff2 = 0;
	    $idA = $ortoA[$i];
	    $idB = $ortoB[$i];
	    $diff1 = $bestscoreAC[$idA] - $ortoS[$i];
	    $diff2 = $bestscoreBC[$idB] - $ortoS[$i];
	    if ($diff1 > $outgroup_cutoff){
				print OUTGR "Ortholog pair $i ($nameA[$idA]-$nameB[$idB]).
   			$nameA[$idA] from $fasta_seq_fileA is closer to $nameC[$besthitAC[$idA]] than to $nameB[$idB]\n";
				print OUTGR "   $ortoS[$i] < $bestscoreAC[$idA] by $diff1\n";
	    }
	    if ($diff2 > $outgroup_cutoff){
				print OUTGR "Ortholog pair $i ($nameA[$idA]-$nameB[$idB]).
   			$nameB[$idB] from $fasta_seq_fileB is closer to $nameC[$besthitBC[$idB]] than to $nameA[$idA]\n";
				print OUTGR "   $ortoS[$i] < $bestscoreBC[$idB] by $diff2\n";
	    }
	    if (($diff1 > $outgroup_cutoff) or ($diff2 > $outgroup_cutoff)){
				++$rejected;
				$del[$i] = 1;
	    }
		}
		print "Number of rejected groups: $rejected (outgroup sequence was closer by more than $outgroup_cutoff bits)\n";
		close OUTGR;
  	move("./$file", "$output_dir/$file");
  } # End of $use_outgroup
################################
# Read inside scores from AA
################################
  $count = 0;
  $max_hit = 0;
	if (-e "$seqSimilarity_outputAA"){
		open AA, "$seqSimilarity_outputAA" or die "Blast output file A->A is missing\n";
	}elsif(-e "$seqSimilarity_outputAA_alternative"){
		open AA, "$seqSimilarity_outputAA_alternative" or die "Blast output file A->A is missing\n";
	}else{
		die "Blast output file A->A is missing\n";
	}
  while (<AA>) {
		chomp;                  # strip newline
		@Fld = split(/\s+/);    # Get query and match names
		if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
				print STDERR "AA ok\n";
	    }
	    next;
		}
		$q = $Fld[0];
		$m = $Fld[1];
		$score = $Fld[2];
		next unless (vec($is_ortologA,$idA{$q},1)); # Removed this line, since we need the scores for all in-paralogs for the consistency check
		next if (!overlap_test(@Fld));
		next if ($score < $score_cutoff);
		if(!$count || $q ne $oldq){ # New query
		    $max_hit = $hit if ($hit > $max_hit);
		    $hit = 0;
		    $oldq = $q;
		}
		++$hit;

		$scoreAA{"$idA{$q}:$idA{$m}"}  = int($score + 0.5);
		$hitAA[$idA{$q}][$hit] = int($idA{$m});
		$hitnAA[$idA{$q}] = $hit;
  	++$count;
  }
  close AA;
  if ($stats){
		print OUTPUT "$count $fasta_seq_fileA_withoutCounter-$fasta_seq_fileA_withoutCounter matches\n";
  }

################################
# Read inside scores from BB
################################
  $count = 0;
	if (-e "$seqSimilarity_outputBB"){
		open BB, "$seqSimilarity_outputBB" or die "Blast output file B->B is missing\n";
	}elsif(-e "$seqSimilarity_outputBB_alternative"){
		open BB, "$seqSimilarity_outputBB_alternative" or die "Blast output file B->B is missing\n";
	}else{
		die "Blast output file B->B is missing\n";
	}
  while (<BB>) {
		chomp;                  # strip newline
		@Fld = split(/\s+/);    # Get query and match names
		if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
				print STDERR "BB ok\n";
	    }
	    next;
		}
		$q = $Fld[0];
		$m = $Fld[1];
		$score = $Fld[2];
		next unless (vec($is_ortologB,$idB{$q},1)); # Removed this line, since we need the scores for all in-paralogs for the consistency check
		next if (!overlap_test(@Fld));
		next if ($score < $score_cutoff);
		if(!$count || $q ne $oldq){ # New query
		    $max_hit = $hit if ($hit > $max_hit);
		    $oldq = $q;
		    $hit = 0;
		}
		++$hit;
		$scoreBB{"$idB{$q}:$idB{$m}"} = int($score + 0.5);
		$hitBB[$idB{$q}][$hit] = int($idB{$m});
		$hitnBB[$idB{$q}] = $hit;
	  ++$count;
  }
  close BB;
  if ($stats){
		print OUTPUT "$count $fasta_seq_fileB_withoutCounter-$fasta_seq_fileB_withoutCounter matches\n";
  }
  if ($show_times){
		($user_time,,,) = times;
		printf ("Reading paralogous hits took %.2f seconds\n", ($user_time - $prev_time));
		$prev_time = $user_time;
  }
  print "Maximum number of hits per sequence was $max_hit\n" if ($debug);

#####################################################
# Find paralogs:
#####################################################
  for $i(1..$o){
		$merge[$i] = 0;
		next if($del[$i]); # If outgroup species was closer to one of the seed orthologs
		$idA = $ortoA[$i];
		$idB = $ortoB[$i];
		local @membersA = ();
		local @membersB = ();
		undef $is_paralogA[$i];
		undef $is_paralogB[$i];
	  print "$i: Ortholog pair $nameA[$idA] and $nameB[$idB]. $hitnAA[$idA] hits for A and $hitnBB[$idB] hits for B\n"  if ($debug);
		# Check if current ortholog is already clustered:
		for $j(1..($i-1)){
	    # Overlap type 1: Both orthologs already clustered here -> merge
	    if ((vec($is_paralogA[$j],$idA,1)) and (vec($is_paralogB[$j],$idB,1))){
				$merge[$i] = $j;
				print "Merge CASE 1: group $i ($nameB[$idB]-$nameA[$idA]) and $j ($nameB[$ortoB[$j]]-$nameA[$ortoA[$j]])\n" if ($debug);
				last;
	    }
	    # Overlap type 2: 2 competing ortholog pairs -> merge
	    elsif (($ortoS[$j] - $ortoS[$i] <= $grey_zone) and (($ortoA[$j] == $ortoA[$i]) or ($ortoB[$j] == $ortoB[$i]))){ # The last condition is false if the previous cluster has been already deleted
				$merge[$i] = $j;
				print "Merge CASE 2: group $i ($nameA[$ortoA[$i]]-$nameB[$ortoB[$i]]) and $j ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
				last;
	    }
	    # Overlap type 3: DELETE One of the orthologs belongs to some much stronger cluster -> delete
	    elsif (((vec($is_paralogA[$j],$idA,1)) or (vec($is_paralogB[$j],$idB,1))) and ($ortoS[$j] - $ortoS[$i] > $score_cutoff)){
				print "Delete CASE 3: Cluster $i -> $j, score $ortoS[$i] -> $ortoS[$j], ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
				$merge[$i]= -1; # Means - do not add sequences to this cluster
				$paralogsA[$i] = "";
				$paralogsB[$i] = "";
				last;
	    }
	    # Overlap type 4: One of the orthologs is close to the center of other cluster
	    elsif (((vec($is_paralogA[$j],$idA,1)) and ($confPA[$idA] > $group_overlap_cutoff)) or ((vec($is_paralogB[$j],$idB,1)) and ($confPB[$idB] > $group_overlap_cutoff))){
				print "Merge CASE 4: Cluster $i -> $j, score $ortoS[$i] -> $ortoS[$j], ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
				$merge[$i] = $j;
				last;
	    }
	    # Overlap type 5:
	    # All clusters that were overlapping, but not catched by previous "if" statements will be DIVIDED!
		}
		next if ($merge[$i] < 0); # This cluster should be deleted
##### Check for paralogs in A
		$N = $hitnAA[$idA];
		for $j(1..$N){
	    $hitID = $hitAA[$idA][$j]; # hit of idA
	    # Decide whether this hit is inside the paralog circle:
	    if ( ($idA == $hitID) or ($scoreAA{"$idA:$hitID"} >= $bestscoreAB[$idA]) and ($scoreAA{"$idA:$hitID"} >= $bestscoreAB[$hitID])){
				if ($debug == 2){
		    	print "   Paralog candidates: ";
		    	printf ("%-20s: %-20s", $nameA[$idA], $nameA[$hitID]);
		    	print "\t$scoreAA{\"$idA:$hitID\"} : $bestscoreAB[$idA] : $bestscoreAB[$hitID]\n";
				}
				$paralogs = 1;
				if ($scoreAA{"$idA:$idA"} == $ortoS[$i]){
		    	if ($scoreAA{"$idA:$hitID"} == $scoreAA{"$idA:$idA"}){
						$conf_here = 1.0; # In the center
		    	}
		    else{
					$conf_here = 0.0; # On the border
		    }
			} else {
		  	$conf_here = ($scoreAA{"$idA:$hitID"} - $ortoS[$i]) / ($scoreAA{"$idA:$idA"} - $ortoS[$i]);
			}
		# Check if this paralog candidate is already clustered in other clusters
			for $k(1..($i-1)){
		    if (vec($is_paralogA[$k],$hitID,1)){ # Yes, found in cluster $k
					if($debug == 2){
				    print "      $nameA[$hitID] is already in cluster $k, together with:";
				    print " $nameA[$ortoA[$k]] and $nameB[$ortoB[$k]] ";
				    print "($scoreAA{\"$ortoA[$k]:$hitID\"})";
					}
					if (($confPA[$hitID] >= $conf_here) and
				    ($j != 1)){ # The seed ortholog CAN NOT remain there
				    print " and remains there.\n" if ($debug == 2);
				    $paralogs = 0; # No action
					}
					else { # Ortholog of THIS cluster is closer than ortholog of competing cluster $k
					    print " and should be moved here!\n" if ($debug == 2); # Remove from other cluster, add to this cluster
					    @membersAK = split(/ /, $paralogsA[$k]); # This array contains IDs
					    $paralogsA[$k] = "";# Remove all paralogs from cluster $k
							@tmp = ();
					    for $m(@membersAK){
								push(@tmp,$m) if ($m != $hitID); # Put other members back
					    }
					    $paralogsA[$k] = join(' ',@tmp);
					    undef $is_paralogA[$k]; # Create index that we can check later
					    grep (vec($is_paralogA[$k],$_,1) = 1, @tmp);
					}
					last;
		    }
			}
			next if (! $paralogs); # Skip that paralog - it is already in cluster $k
			push (@membersA,$hitID); # Add this hit to paralogs of A
	  }
	}
	# Calculate confidence values now:
	@tmp = ();
	for $idP (@membersA){ # For each paralog calculate conf value
	    if($scoreAA{"$idA:$idA"} == $ortoS[$i]){
				if ($scoreAA{"$idA:$idP"} == $scoreAA{"$idA:$idA"}){
		    	$confPA[$idP] = 1.00;
				}
				else{
		    	$confPA[$idP] = 0.00;
				}
	    }
	    else{
				$confPA[$idP] = ($scoreAA{"$idA:$idP"} - $ortoS[$i]) / ($scoreAA{"$idA:$idA"} - $ortoS[$i]);
	    }
	    push (@tmp, $idP) if ($confPA[$idP] >= $conf_cutoff); # If one wishes to use only significant paralogs
	}

  @membersA = @tmp;
	########### Merge if necessary:
	if ($merge[$i] > 0){ # Merge existing cluster with overlapping cluster
	    @tmp = split(/ /,$paralogsA[$merge[$i]]);
	    for $m (@membersA){
				push (@tmp, $m)  unless (vec($is_paralogA[$merge[$i]],$m,1));
	    }
	    $paralogsA[$merge[$i]] = join(' ',@tmp);
	    undef $is_paralogA[$merge[$i]];
	    grep (vec($is_paralogA[$merge[$i]],$_,1) = 1, @tmp); # Refresh index of paralog array
	}
	######### Typical new cluster:
	else{  # Create a new cluster
	    $paralogsA[$i] = join(' ',@membersA);
	    undef $is_paralogA; # Create index that we can check later
	    grep (vec($is_paralogA[$i],$_,1) = 1, @membersA);
	}

##### The same procedure for species B:
	$N = $hitnBB[$idB];
	for $j(1..$N){
	    $hitID = $hitBB[$idB][$j];
	    if ( ($idB == $hitID) or ($scoreBB{"$idB:$hitID"} >= $bestscoreBA[$idB]) and ($scoreBB{"$idB:$hitID"} >= $bestscoreBA[$hitID])){
				if ($debug == 2){
			    print "   Paralog candidates: ";
			    printf ("%-20s: %-20s", $nameB[$idB], $nameB[$hitID]);
			    print "\t$scoreBB{\"$idB:$hitID\"} : ";
			    print "$bestscoreBA[$idB] : $bestscoreBA[$hitID]\n";
				}
				$paralogs = 1;
				if ($scoreBB{"$idB:$idB"} == $ortoS[$i]){
				    if ($scoreBB{"$idB:$hitID"} == $scoreBB{"$idB:$idB"}){
							$conf_here = 1.0;
				    }
				    else{
							$conf_here = 0.0;
				    }
				}
				else{
				    $conf_here = ($scoreBB{"$idB:$hitID"} - $ortoS[$i]) / ($scoreBB{"$idB:$idB"} - $ortoS[$i]);
				}

		# Check if this paralog candidate is already clustered in other clusters
		for $k(1..($i-1)){
		    if (vec($is_paralogB[$k],$hitID,1)){ # Yes, found in cluster $k
					if($debug == 2){
				    print "      $nameB[$hitID] is already in cluster $k, together with:";
				    print " $nameB[$ortoB[$k]] and $nameA[$ortoA[$k]] ";
				    print "($scoreBB{\"$ortoB[$k]:$hitID\"})";
					}
					if (($confPB[$hitID] >= $conf_here) and
				    ($j != 1)){ # The seed ortholog CAN NOT remain there
				    print " and remains there.\n" if ($debug == 2);
				    $paralogs = 0; # No action
					}
					else { # Ortholog of THIS cluster is closer than ortholog of competing cluster $k
				    print " and should be moved here!\n" if ($debug == 2); # Remove from other cluster, add to this cluster
				    @membersBK = split(/ /, $paralogsB[$k]); # This array contains names, not IDs
				    $paralogsB[$k] = "";
				    @tmp = ();
			    	for $m(@membersBK){
							push(@tmp,$m) if ($m != $hitID); # Put other members back
			    	}
			    	$paralogsB[$k] = join(' ',@tmp);
			    	undef $is_paralogB[$k]; # Create index that we can check later
			    	grep (vec($is_paralogB[$k],$_,1) = 1, @tmp);
					}
					last; # Don't search in other clusters
		    }
		}
		next if (! $paralogs); # Skip that paralog - it is already in cluster $k
		push (@membersB,$hitID);
	 	}
	}
	# Calculate confidence values now:
	@tmp = ();
	for $idP (@membersB){ # For each paralog calculate conf value
	    if($scoreBB{"$idB:$idB"} == $ortoS[$i]){
				if ($scoreBB{"$idB:$idP"} == $scoreBB{"$idB:$idB"}){
		    	$confPB[$idP] = 1.0;
				}
				else{
		    	$confPB[$idP] = 0.0;
				}
	    }
	    else{
				$confPB[$idP] = ($scoreBB{"$idB:$idP"} - $ortoS[$i]) / ($scoreBB{"$idB:$idB"} - $ortoS[$i]);
	    }
	    push (@tmp, $idP) if ($confPB[$idP] >= $conf_cutoff); # If one wishes to use only significant paralogs
	}

	@membersB = @tmp;
	########### Merge if necessary:
	if ($merge[$i] > 0){ # Merge existing cluster with overlapping cluster
	    @tmp = split(/ /,$paralogsB[$merge[$i]]);
	    for $m (@membersB){
				push (@tmp, $m)  unless (vec($is_paralogB[$merge[$i]],$m,1));
	    }
	    $paralogsB[$merge[$i]] = join(' ',@tmp);
	    undef $is_paralogB[$merge[$i]];
	    grep (vec($is_paralogB[$merge[$i]],$_,1) = 1, @tmp); # Refresh index of paralog array
	}
	######### Typical new cluster:
	else{  # Create a new cluster
	    $paralogsB[$i] = join(' ',@membersB);
	    undef $is_paralogB; # Create index that we can check later
	    grep (vec($is_paralogB[$i],$_,1) = 1, @membersB);
	}
}
if ($show_times){
	($user_time,,,) = times;
	printf ("Finding in-paralogs took %.2f seconds\n", ($user_time - $prev_time));
	$prev_time = $user_time;
}
#####################################################
  clean_up(1);
####################################################
# Find group for orphans. If cluster contains only one member, find where it should go:
for $i (1..$o){
	@membersA = split(/ /, $paralogsA[$i]);
	@membersB = split(/ /, $paralogsB[$i]);
	$na = @membersA;
	$nb = @membersB;
	if (($na == 0) and $nb){
		print "Warning: empty A cluster $i\n";
	  for $m (@membersB){
			$bestscore = 0;
			$bestgroup = 0;
			$bestmatch = 0;
			for $j (1..$o) {
		    next if ($i == $j); # Really need to check against all 100% members of the group.
		    @membersBJ = split(/ /, $paralogsB[$j]);
		    for $k (@membersBJ){
					next if ($confPB[$k] != 1); # For all 100% in-paralogs
					$score = $scoreBB{"$m:$k"};
					if ($score > $bestscore){
			    	$bestscore = $score;
			    	$bestgroup = $j;
			    	$bestmatch = $k;
					}
		    }
			}
			print "Orphan $nameB[$m] goes to group $bestgroup with $nameB[$bestmatch]\n" ;
			@members = split(/ /, $paralogsB[$bestgroup]);
			push (@members, $m);
			$paralogsB[$bestgroup] = join(' ',@members);
			$paralogsB[$i] = "";
			undef $is_paralogB[$bestgroup];
			undef $is_paralogB[$i];
			grep (vec($is_paralogB[$bestgroup],$_,1) = 1, @members); # Refresh index of paralog array
	  }
	}
	if ($na and ($nb == 0)){
	  print "Warning: empty B cluster $i\n";
	  for $m (@membersA){
			$bestscore = 0;
			$bestgroup = 0;
			$bestmatch = 0;
			for $j (1..$o) {
		    next if ($i == $j);
		    @membersAJ = split(/ /, $paralogsA[$j]);
		    for $k (@membersAJ){
					next if ($confPA[$k] != 1); # For all 100% in-paralogs
					$score = $scoreAA{"$m:$k"};
					if ($score > $bestscore){
			    	$bestscore = $score;
			    	$bestgroup = $j;
			    	$bestmatch = $k;
					}
		   }
		}
		print "Orphan $nameA[$m] goes to group $bestgroup with $nameA[$bestmatch]\n";
		@members = split(/ /, $paralogsA[$bestgroup]);
		push (@members, $m);
		$paralogsA[$bestgroup] = join(' ',@members);
		$paralogsA[$i] = "";
		undef $is_paralogA[$bestgroup];
		undef $is_paralogA[$i];
		grep (vec($is_paralogA[$bestgroup],$_,1) = 1, @members); # Refresh index of paralog array
	}
}
}

clean_up(1);
my $groupCount=0;
for $i(1..$o){
	@members = split(/ /, $paralogsA[$i]);
	if (@members){
		$groupCount+=1;
	}
}
###################
$htmlfile = "orthologs." . $fasta_seq_fileA . "-" . $fasta_seq_fileB . ".html";
if ($html){
	open HTML, ">$htmlfile" or warn "Could not write to HTML file $filename\n";
}

if ($stats){
	print OUTPUT "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print OUTPUT "$groupCount groups of orthologs\n";
	print OUTPUT "$totalA in-paralogs from $fasta_seq_fileA_withoutCounter\n";
	print OUTPUT "$totalB in-paralogs from $fasta_seq_fileB_withoutCounter\n";
	print OUTPUT "Grey zone $grey_zone bits\n";
	print OUTPUT "Score cutoff $score_cutoff bits\n";
	print OUTPUT "In-paralogs with confidence less than $conf_cutoff not shown\n";
	print OUTPUT "Sequence overlap cutoff $seq_overlap_cutoff\n";
	print OUTPUT "Group merging cutoff $group_overlap_cutoff\n";
	print OUTPUT "Scoring matrix $matrix\n";
	print OUTPUT "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
}
if ($html){
	print HTML "<pre>\n";
	print HTML "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print HTML "$groupCount groups of orthologs\n";
	print HTML "$totalA in-paralogs from $fasta_seq_fileA_withoutCounter\n";
	print HTML "$totalB in-paralogs from $fasta_seq_fileB_withoutCounter\n";
	print HTML "Grey zone $grey_zone bits\n";
	print HTML "Score cutoff $score_cutoff bits\n";
	print HTML "In-paralogs with confidence less than $conf_cutoff not shown\n";
	print HTML "Sequence overlap cutoff $seq_overlap_cutoff\n";
	print HTML "Group merging cutoff $group_overlap_cutoff\n";
	print HTML "Scoring matrix $matrix\n";
	print HTML "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
}
# ##############################################################################
# Check for alternative orthologs, sort paralogs by confidence and print results
# ##############################################################################
if ($use_bootstrap and $debug){
	open FF, ">BS_vs_bits" or warn "Could not write to file BS_vs_bits\n";
}
my $groupCount=0;
for $i(1..$o){
	@membersA = split(/ /, $paralogsA[$i]);
	my $emptyGroup=0;
	if (@membersA){
		$groupCount+=1;
	}else{
		$emptyGroup=1;
	}
	@membersB = split(/ /, $paralogsB[$i]);
	$message = "";
	$htmlmessage = "";

	$idB = $ortoB[$i];
	$nB = $hitnBA[$idB];
	for $idA(@membersA){
	    next if ($confPA[$idA] != 1.0);
	    $nA = $hitnAB[$idA];
	    $confA[$i] = $ortoS[$i]; # default
	    $bsA[$idA] = 1.0;
			$seedScoreA[$idA]=1.0;
			$nextScoreA[$idA]=0;
			$diffNextScoreA[$idA]=0;
	    ##############
	    for $j(1..$nB){
				$idH = $hitBA[$idB][$j];
				################ Some checks for alternative orthologs:
				# 1. Don't consider sequences that are already in this cluster
				next if (vec($is_paralogA[$i],$idH,1));
				next if ($confPA[$idH] > 0); # If $conf_cutoff > 0 idH might be incide circle, but not paralog

				# 2. Check if candidate for alternative ortholog is already clustered in stronger clusters
				$in_other_cluster = 0;
				for $k(1..($i-1)){ # Check if current ortholog is already clustered
		    	if (vec($is_paralogA[$k],$idH,1)){
						$in_other_cluster = $k;
						last;
		    	}
				}

				# 3. The best hit of candidate ortholog should be ortoA or at least to belong into this cluster
				@besthits = split (/ /,$besthitAB[$idH]);
				$this_family = 0;
				for $bh (@besthits){
		    	$this_family = 1 if ($idB == $bh);
				}

				################# Done with checks - if sequence passed, then it could be an alternative ortholog
				$confA[$i] = $ortoS[$i] - $scoreBA{"$idB:$idH"};
				if ($use_bootstrap){
		    	if ($confA[$i] < $ortoS[$i]){ # Danger zone - check for bootstrap
						$bsA[$idA] = bootstrap($fasta_seq_fileB,$idB,$idA,$idH);
		    	}
		    	else {
						$bsA[$idA] = 1.0;
		    	}
				}
				$seedScoreA[$idA] = ( (1.32250111*$confA[$i]) - (0.0110025827*($confA[$i]**2)) + (0.0000294147407*($confA[$i]**3)) + 47.884217)/100;
				if($seedScoreA[$idA] >1.0){
					$seedScoreA[$idA] =1.0
				}
				$nextScoreA[$idA]=$scoreBA{"$idB:$idH"};
				$diffNextScoreA[$idA]=$confA[$i];
				last;
	    }
	    $message .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.", $nameA[$idA], 100*$bsA[$idA]);
	    $message .= sprintf(" Alternative seed ortholog is %s (%d bits away from this cluster)", $nameA[$idH], $confA[$i]) if ($bsA[$idA] < 0.75);
	    $message .= sprintf("\n");
	    if ($html){
				if ($bsA[$idA] < 0.75){
		    	$htmlmessage .= sprintf("<font color=\"red\">");
				}
				elsif ($bsA[$idA] < 0.95){
		    	$htmlmessage .= sprintf("<font color=\"\#FFCC00\">");
				}
				else {
		    	$htmlmessage .= sprintf("<font color=\"green\">");
				}
				$htmlmessage .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.\n", $nameA[$idA], 100*$bsA[$idA]);
				$htmlmessage .= sprintf("Alternative seed ortholog is %s (%d bits away from this cluster)\n", $nameA[$idH], $confA[$i]) if ($bsA[$idA] < 0.75);
				$htmlmessage .= sprintf("</font>");
	   }
	   printf (FF "%s\t%d\t%d\n", $nameA[$idA], $confA[$i], 100*$bsA[$idA]) if ($use_bootstrap and $debug);
	}
	########
	$idA = $ortoA[$i];
	$nA = $hitnAB[$idA];
	for $idB(@membersB){
	    next if ($confPB[$idB] != 1.0);
	    $nB = $hitnBA[$idB];
	    $confB[$i] = $ortoS[$i]; # default
	    $bsB[$idB] = 1.0;
			$seedScoreB[$idB]=1.0;
			$nextScoreB[$idB]=0;
			$diffNextScoreB[$idB]=0;
	    for $j(1..$nA){ # For all AB hits of given ortholog
				$idH = $hitAB[$idA][$j];
				# ############### Some checks for alternative orthologs:
				# 1. Don't consider sequences that are already in this cluster
				next if (vec($is_paralogB[$i],$idH,1));
				next if ($confPB[$idH] > 0); # If $conf_cutoff > 0 idH might be incide circle, but not paralog

				# 2. Check if candidate for alternative ortholog is already clustered in stronger clusters
				$in_other_cluster = 0;
				for $k(1..($i-1)){
		    	if (vec($is_paralogB[$k],$idH,1)){
						$in_other_cluster = $k;
						last; # out from this cycle
		    	}
				}

				# 3. The best hit of candidate ortholog should be ortoA
				@besthits = split (/ /,$besthitBA[$idH]);
				$this_family = 0;
				for $bh (@besthits){
		    	$this_family = 1 if ($idA == $bh);
				}

				# ################ Done with checks - if sequence passed, then it could be an alternative ortholog
				$confB[$i] = $ortoS[$i] - $scoreAB{"$idA:$idH"};
				if ($use_bootstrap){
		    	if ($confB[$i] < $ortoS[$i]){
						$bsB[$idB] = bootstrap($fasta_seq_fileA,$idA,$idB,$idH);
		    	}
		    	else {
						$bsB[$idB] = 1.0;
		    	}
				}
				$seedScoreB[$idB] = (40.482*($confB[$i]**0.201))/100;
				if($seedScoreB[$idB]>1.0){
					$seedScoreB[$idB]=1.0
				}
				$nextScoreB[$idB]=$scoreAB{"$idA:$idH"};
				$diffNextScoreB[$idB]=$confB[$i];
				last;
	    }
	    $message .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.", $nameB[$idB], 100*$bsB[$idB]);
	    $message .= sprintf(" Alternative seed ortholog is %s (%d bits away from this cluster)", $nameB[$idH],$confB[$i]) if ($bsB[$idB] < 0.75);
	    $message .= sprintf("\n");
	    if ($html){
				if ($bsB[$idB] < 0.75){
		    	$htmlmessage .= sprintf("<font color=\"red\">");
				}
				elsif ($bsB[$idB] < 0.95){
		    	$htmlmessage .= sprintf("<font color=\"\#FFCC00\">");
				}
				else {
		    	$htmlmessage .= sprintf("<font color=\"green\">");
				}
				$htmlmessage .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.\n", $nameB[$idB], 100*$bsB[$idB]);
				$htmlmessage .= sprintf("Alternative seed ortholog is %s (%d bits away from this cluster)\n",  $nameB[$idH],$confB[$i]) if ($bsB[$idB] < 0.75);
				$htmlmessage .= sprintf("</font>");
	    }
	    printf (FF "%s\t%d\t%d\n", $nameB[$idB], $confB[$i], 100*$bsB[$idB]) if ($use_bootstrap and $debug);
		}
		close FF;
	########### Print header ###############
		if ($stats){
			if($emptyGroup==0){
				print OUTPUT "___________________________________________________________________________________\n";
		 		print OUTPUT "Group of orthologs #" . $groupCount .". Best score $ortoS[$i] bits\n";
	 			print OUTPUT "Score difference with first non-orthologous sequence - ";
		 		printf (OUTPUT "%s:%d   %s:%d\n", $fasta_seq_fileA_withoutCounter,$confA[$i],$fasta_seq_fileB_withoutCounter,$confB[$i]);
			}
		}

		if ($html){
			if($emptyGroup==0){
				print HTML "</pre>\n";
		 		print HTML "<hr WIDTH=\"100%\">";
		 		print HTML "<h3>";
		 		print HTML "Group of orthologs #" . $groupCount .". Best score $ortoS[$i] bits<br>\n";
	 			print HTML "Score difference with first non-orthologous sequence - ";
		 		printf (HTML "%s:%d   %s:%d</h3><pre>\n", $fasta_seq_fileA_withoutCounter,$confA[$i],$fasta_seq_fileB_withoutCounter,$confB[$i]);
			}
		}

	########### Sort and print members of A ############
  @membersAOrderedByScore=();
  %paralogsScoreHash=();
  for $m (@membersA){
    $paralogsScoreHash{$m}=$confPA[$m];
  }
  foreach my $p (sort { $paralogsScoreHash{$b} <=> $paralogsScoreHash{$a} } keys %paralogsScoreHash) {
    push (@membersAOrderedByScore, $p);
  }
  @membersA=@membersAOrderedByScore;
	$paralogsA[$i] = join(' ',@membersA); # Put them back together

	# Sort membersB inside the cluster by confidence:
  @membersBOrderedByScore=();
  %paralogsScoreHash=();
  for $m (@membersB){
    $paralogsScoreHash{$m}=$confPB[$m];
  }
  foreach my $p (sort { $paralogsScoreHash{$b} <=> $paralogsScoreHash{$a} } keys %paralogsScoreHash) {
    push (@membersBOrderedByScore, $p);
  }
  @membersB=@membersBOrderedByScore;
	$paralogsB[$i] = join(' ',@membersB); # Put them back together

  $nA = @membersA;
	$nB = @membersB;
	$nMAX = ($nA > $nB) ? $nA : $nB;
  # Print to text file and to HTML file
	for $m (0..($nMAX-1)){
	    if ($m < $nA){
				if ($stats){
		    	printf (OUTPUT "%-20s\t%.2f%%\t\t", $nameA[$membersA[$m]], (100*$confPA[$membersA[$m]]));
				}
				if ($html){
		    	print HTML "<B>" if ($confPA[$membersA[$m]] == 1);
		    	printf (HTML "%-20s\t%.2f%%\t\t", $nameA[$membersA[$m]], (100*$confPA[$membersA[$m]]));
		    	print HTML "</B>" if ($confPA[$membersA[$m]] == 1);
				}
	    }
	    else {
				printf (OUTPUT "%-20s\t%-7s\t\t", "                    ", "       ");
				printf (HTML "%-20s\t%-7s\t\t", "                    ", "       ") if ($html);
	    }
	    if ($m < $nB){
				if ($stats){
		    	printf (OUTPUT "%-20s\t%.2f%%\n", $nameB[$membersB[$m]], (100*$confPB[$membersB[$m]]));
				}
				if ($html){
		    	print HTML "<B>" if ($confPB[$membersB[$m]] == 1);
		    	printf (HTML "%-20s\t%.2f%%", $nameB[$membersB[$m]], (100*$confPB[$membersB[$m]]));
		    	print HTML "</B>" if ($confPB[$membersB[$m]] == 1);
		    	print HTML "\n";
				}
	    }
	    else {
				printf (OUTPUT "%-20s\t%-7s\n", "                    ", "       ") if($stats);
				print HTML "\n" if ($html);
	    }
		}
		print OUTPUT $message if ($use_bootstrap and $stats);
		print HTML "$htmlmessage" if ($use_bootstrap and $html);
  }
  if ($stats) {
  	close OUTPUT;
    move("./$statsfile", "$output_dir/$statsfile");
    print "Output saved to file $output_dir/$statsfile\n";
  }
  if ($html){
    close HTML;
    move("./$htmlfile", "$output_dir/$htmlfile");
    print "HTML output saved to $output_dir/$htmlfile\n";
  }
  if ($table){
		$filename = "table." . $fasta_seq_fileA . "-" . $fasta_seq_fileB;
		open F, ">$filename" or die;
		print F "OrtoID\tScore\tOrtoA\tOrtoB\n";
		my $groupCount=0;
		for $i(1..$o){
				@members = split(/ /, $paralogsA[$i]);
				if (@members){
					$groupCount+=1;
					print F "$groupCount\t$ortoS[$i]\t";
				}
			  for $m (@members){
					$m =~ s/://g;
					printf (F "%s %.3f ", $nameA[$m], $confPA[$m]);
			  }
				if (@members){
			    print F "\t";
				}
			  @members = split(/ /, $paralogsB[$i]);
			  for $m (@members){
					$m =~ s/://g;
					printf (F "%s %.3f ", $nameB[$m], $confPB[$m]);
			  }
				if (@members){
			    print F "\n";
				}
			}
			close F;
		  move("./$filename", "$output_dir/$filename");
			print "Table output saved to $output_dir/$filename\n";
    }
    if ($mysql_table){
			$filename2 = "SQLtable." . $fasta_seq_fileA . "-" . $fasta_seq_fileB;
			open F2, ">$filename2" or die;
			my $groupCount=0;
			for $i(1..$o){
			    @membersA = split(/ /, $paralogsA[$i]);
					if (@membersA){
						$groupCount+=1;
					}
			    for $m (@membersA){
						if ($use_bootstrap && $bsA[$m]) {
							if ($seedscore && $seedScoreA[$m]){
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\t%.3f\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m],100*$bsA[$m], $seedScoreA[$m]);
								#printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\t%.3f\t%d\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m],100*$bsA[$m], $seedScoreA[$m],$diffNextScoreA[$m]);
							}else{
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m], 100*$bsA[$m]);
							}
						} else {
							if ($seedscore && $seedScoreA[$m]){
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%.3f\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m], $seedScoreA[$m]);
								#printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%.3f\t%d\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m], $seedScoreA[$m],$diffNextScoreA[$m]);
							}else{
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\n", $groupCount, $ortoS[$i], $fasta_seq_fileA_withoutCounter, $confPA[$m], $nameA[$m]);
							}
						}
			    }
			    @membersB = split(/ /, $paralogsB[$i]);
			    for $m (@membersB){
						if ($use_bootstrap && $bsB[$m]) {
							if ($seedscore && $seedScoreB[$m]){
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\t%.3f\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m], 100*$bsB[$m], $seedScoreB[$m]);
								#printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\t%.3f\t%d\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m], 100*$bsB[$m], $seedScoreB[$m],$diffNextScoreB[$m]);
							}else{
								printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%d%\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m], 100*$bsB[$m]);
							}
						}else {
								if ($seedscore && $seedScoreB[$m]){
									printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%.3f\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m],  $seedScoreB[$m]);
									#printf (F2 "%d\t%d\t%s\t%.3f\t%s\t%.3f\t%d\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m],  $seedScoreB[$m],$diffNextScoreB[$m]);
								}else{
									printf (F2 "%d\t%d\t%s\t%.3f\t%s\n", $groupCount, $ortoS[$i], $fasta_seq_fileB_withoutCounter, $confPB[$m], $nameB[$m]);
								}
						}
			    }
			}
			close (F2);
		  move("./$filename2", "$output_dir/$filename2");
			print "mysql output saved to $output_dir/$filename2\n";
    }
    if ($show_times){
      ($user_time,,,) = times;
      printf ("Finding bootstrap values and printing took %.2f seconds\n", ($user_time - $prev_time));
    }
    if (0 && $run_blast) {
      unlink "formatdb.log";
      unlink "$fasta_seq_fileA.phr", "$fasta_seq_fileA.pin", "$fasta_seq_fileA.psq";
      unlink "$fasta_seq_fileB.phr", "$fasta_seq_fileB.pin", "$fasta_seq_fileB.psq";
      unlink "$fasta_seq_fileC.phr", "$fasta_seq_fileC.pin", "$fasta_seq_fileC.psq" if ($use_outgroup);
    }
  }

##############################################################
# Functions:
##############################################################
sub clean_up { # Sort members within cluster and clusters by size
############################################################################################### Modification by Isabella 3
	my @Fld = @_;
	# Sort on index arrays with perl's built in sort instead of using bubble sort.
  my $var = $Fld[0];
  $totalA = $totalB = 0;
  # First pass: count members within each cluster
  foreach $i (1..$o) {
		@membersA = split(/ /, $paralogsA[$i]);
		$clusnA[$i] = @membersA; # Number of members in this cluster
		$totalA += $clusnA[$i];
		$paralogsA[$i] = join(' ',@membersA);
		@membersB = split(/ /, $paralogsB[$i]);
		$clusnB[$i] = @membersB; # Number of members in this cluster
		$totalB += $clusnB[$i];
		$paralogsB[$i] = join(' ',@membersB);

		$clusn[$i] =  $clusnB[$i] + $clusnA[$i]; # Number of members in given group
  }

  # Create an array used to store the position each element shall have in the final array
  # The elements are initialized with the position numbers
  my @position_index_array = (1..$o);

  # Sort the position list according to cluster size
  my @cluster_sorted_position_list = sort { $clusn[$b] <=> $clusn[$a]} @position_index_array;

  # Create new arrays for the sorted information
  my @new_paralogsA;
  my @new_paralogsB;
  my @new_is_paralogA;
  my @new_is_paralogB;
  my @new_clusn;
  my @new_ortoS;
  my @new_ortoA;
  my @new_ortoB;
  # Add the information to the new arrays in the orer specifeid by the index array
  for (my $index_in_list = 0; $index_in_list < scalar @cluster_sorted_position_list; $index_in_list++) {
		my $old_index = $cluster_sorted_position_list[$index_in_list];
		if (!$clusn[$old_index]) {
			$o = (scalar @new_ortoS) - 1;
			last;
		}
		$new_paralogsA[$index_in_list + 1] = $paralogsA[$old_index];
    $new_paralogsB[$index_in_list + 1] = $paralogsB[$old_index];
    $new_is_paralogA[$index_in_list + 1] = $is_paralogA[$old_index];
    $new_is_paralogB[$index_in_list + 1] = $is_paralogB[$old_index];
   	$new_clusn[$index_in_list + 1] = $clusn[$old_index];
		$new_ortoA[$index_in_list + 1] = $ortoA[$old_index];
		$new_ortoB[$index_in_list + 1] = $ortoB[$old_index];
		$new_ortoS[$index_in_list + 1] = $ortoS[$old_index];
	}

  @paralogsA = @new_paralogsA;
  @paralogsB = @new_paralogsB;
  @is_paralogA = @new_is_paralogA;
  @is_paralogB = @new_is_paralogB;
  @clusn = @new_clusn;
  @ortoS = @new_ortoS;
  @ortoA = @new_ortoA;
  @ortoB = @new_ortoB;

  # Create an array used to store the position each element shall have in the final array
  # The elements are initialized with the position numbers
  @position_index_array = (1..$o);

  # Sort the position list according to score
  @score_sorted_position_list = sort { $ortoS[$b] <=> $ortoS[$a] } @position_index_array;

  # Create new arrays for the sorted information
  my @new_paralogsA2 = ();
  my @new_paralogsB2 = ();
  my @new_is_paralogA2 = ();
  my @new_is_paralogB2 = ();
  my @new_clusn2 = ();
  my @new_ortoS2 = ();
  my @new_ortoA2 = ();
  my @new_ortoB2 = ();

  # Add the information to the new arrays in the orer specifeid by the index array
  for (my $index_in_list = 0; $index_in_list < scalar @score_sorted_position_list; $index_in_list++) {
		my $old_index = $score_sorted_position_list[$index_in_list];
		$new_paralogsA2[$index_in_list + 1] = $paralogsA[$old_index];
    $new_paralogsB2[$index_in_list + 1] = $paralogsB[$old_index];
    $new_is_paralogA2[$index_in_list + 1] = $is_paralogA[$old_index];
    $new_is_paralogB2[$index_in_list + 1] = $is_paralogB[$old_index];
   	$new_clusn2[$index_in_list + 1] = $clusn[$old_index];
		$new_ortoA2[$index_in_list + 1] = $ortoA[$old_index];
		$new_ortoB2[$index_in_list + 1] = $ortoB[$old_index];
		$new_ortoS2[$index_in_list + 1] = $ortoS[$old_index];
  }

  @paralogsA = @new_paralogsA2;
  @paralogsB = @new_paralogsB2;
  @is_paralogA = @new_is_paralogA2;
  @is_paralogB = @new_is_paralogB2;
  @clusn = @new_clusn2;
  @ortoS = @new_ortoS2;
  @ortoA = @new_ortoA2;
  @ortoB = @new_ortoB2;

#################################################################################### End modification by Isabella 3
}
sub bootstrap{
  my @Fld = @_;
	my $species = $Fld[0];
  my $seq_id1 = $Fld[1]; # Query ID from $species
  my $seq_id2 = $Fld[2]; # Best hit ID from other species
  my $seq_id3 = $Fld[3]; # Second best hit
  # Retrieve sequence 1 from $species and sequence 2 from opposite species
  my $significance = 0.0;

  if ($species eq $fasta_seq_fileA){
		$file1 = $fasta_seq_fileA;
		$file2 = $fasta_seq_fileB;
  } elsif ($species eq $fasta_seq_fileB){
		$file1 = $fasta_seq_fileB;
		$file2 = $fasta_seq_fileA;
  } else {
		print "Bootstrap values for ortholog groups are not calculated\n";
		return 0;
  }

  open A, $file1 or die;
  $id = 0;
  $print_this_seq = 0;
  $seq1 = "";
  $seq2 = "";

 # $query_file = $seq_id1 . ".faq";
 $query_file = "$seq_id1-$file1.faq";
  open Q, ">$query_file" or die;

  while (<A>){
		if(/^\>/){
	    ++$id;
	    $print_this_seq = ($id == $seq_id1)?1:0;
		}
		print Q if ($print_this_seq);
  }
  close A;
  close Q;
  ###
  open B, $file2 or die;
 # $db_file = $seq_id2 . ".fas";
  $db_file = "$seq_id2-$file2.fas";
  open DB, ">$db_file" or die;
  $id = 0;
  $print_this_seq = 0;
my $parsedFile="$seq_id2-$file2.faa";

  while (<B>){
		if(/^\>/){
	    ++$id;
	    $print_this_seq = (($id == $seq_id2) or ($id == $seq_id3))?1:0;
		}
		print DB if ($print_this_seq);
  }
  close B;
  close DB;
	if ($run_blast){
		system "$formatdb -i $db_file";
		# Use soft masking in 1-pass mode for simplicity.
		system "$blastall -F\"m S\" -i $query_file -z 5000000 -d $db_file -p blastp -M $matrix -m7 | ./$blastParser 0 -a > $parsedFile";
	}elsif ($run_diamond){
		require "$diamondParser";
		system ("$diamond makedb --in $db_file -d $db_file.db --quiet");
		if ($cores>0){
			system ("$diamond blastp -d $db_file.db.dmnd -q $query_file -o diamond$query_file$db_file -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send qseq sseq btop --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0 -p $cores");
		}else{
			system ("$diamond blastp -d $db_file.db.dmnd -q $query_file -o diamond$query_file$db_file -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send qseq sseq btop --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0");
		}
		diamondParser::parse("diamond$query_file$db_file", $score_cutoff, $parsedFile, 0, 1);
    unlink "$db_file.db.dmnd", "diamond$query_file$db_file";
	}
  if (-s($parsedFile) and !-z($parsedFile)){
	#	print $parsedFile;
		system("java -jar $seqstat -m $matrix -n 1000 -i $parsedFile > $seq_id2-$file2.bs"); # Can handle U, u
		if (-s("$seq_id2-$file2.bs")){
	    open BS, "$seq_id2-$file2.bs" or die "pac failed\n";
	    $_ = <BS>;
	    ($dummy1,$dummy2,$dummy3,$dummy4,$significance) = split(/\s+/);
	    close BS;
		}
		else{
	    print STDERR "pac failed\n"; # if ($debug);
	    $significance = -0.01;
		}
  } else{
		print STDERR "blast2faa for $query_file / $db_file failed\n"; # if ($debug);
		$significance = 0.0;
  }

  unlink $db_file, $parsedFile, "$seq_id2-$file2.bs", $query_file;
  unlink "formatdb.log", "$db_file.psq", "$db_file.pin", "$db_file.phr";


  return $significance;
}

sub overlap_test{
  my @Fld = @_;
	# Filter out fragmentary hits by:
	# Ignore hit if aggregate matching area covers less than $seq_overlap_cutoff of sequence.
	# Ignore hit if local matching segments cover less than $segment_coverage_cutoff of sequence.
	# $Fld[3] and $Fld[4] are query and subject lengths.
	# $Fld[5] and $Fld[6] are lengths of the aggregate matching region on query and subject. (From start of first matching segment to end of last matching segment).
	# $Fld[7] and $Fld[8] are local matching length on query and subject (Sum of all segments length's on query).

	$retval = 1;
	if ($Fld[5] < ($seq_overlap_cutoff * $Fld[3])) {$retval = 0};
	if ($Fld[7] < ($segment_coverage_cutoff * $Fld[3])) {$retval = 0};
	if ($Fld[6] < ($seq_overlap_cutoff * $Fld[4])) {$retval = 0};
	if ($Fld[8] < ($segment_coverage_cutoff * $Fld[4])) {$retval = 0};
	return $retval;
}

sub do_blast {
  my @Fld = @_;
  system("touch $Fld[4]");  # So that a parallel job won't start this pair run also.
                            # IMPORTANT note: all parallel jobs must be finised before
                            # starting the clustering phase since blasting may not be finished.
                            # Hence, set run_inparanoid=0 when running parallel blasts.
                            # When all blast jobs are finished, the clustering can be run in parallel.
  print STDERR "Formatting BLAST databases for $Fld[0] and $Fld[1]\n";
  system ("$formatdb -i $Fld[0]") if !-e $Fld[0].".pin";
  system ("$formatdb -i $Fld[1]") if !-e $Fld[1].".pin";

  if ($two_passes) {
    do_blast_2pass(@_);
  } else {
    do_blast_1pass(@_);
  }
  unlink "$Fld[1].psq", "$Fld[0].psq", "$Fld[0].pin", "$Fld[1].pin", "$Fld[0].phr", "$Fld[1].phr", "formatdb.log";
}

sub do_blast_1pass {
  my @Fld = @_;
  # Use soft masking (low complexity masking by SEG in search phase, not in alignment phase).
  print STDERR "Starting BLAST searches...\n";
  system ("$blastall -F\"m S\" -i $Fld[0] -d $Fld[1] -p blastp -v $Fld[3] -b $Fld[3] -M $matrix -z 5000000 -m7 | ./$blastParser $score_cutoff > $Fld[4]");
}

sub do_blast_2pass {
	my @Fld = @_;
	# assume the script has already formatted the database
	# we will now do 2-pass approach
	# load sequences
	%sequencesA = ();
	%sequencesB = ();

	open (FHA, $Fld [0]);
	while (<FHA>) {

		$aLine = $_;
		chomp ($aLine);
		$seq = "";
		if ($aLine =~ />/) {
			$header = $aLine;
			$sequencesA {$header} = "";
		}
		else {
			$sequencesA {$header} = $sequencesA {$header}.$aLine;
		}
	}
	close (FHA);

	open (FHB, $Fld [1]);
	while (<FHB>) {
		$aLine = $_;
		chomp ($aLine);
		$seq = "";
		if ($aLine =~ />/) {
			$header = $aLine;
			$sequencesB {$header} = "";
		}
		else {
			$sequencesB {$header} = $sequencesB {$header}.$aLine;
		}
	}
	close (FHB);

	# Do first pass with compositional adjustment on and soft masking.
	# This efficiently removes low complexity matches but truncates alignments,
	# making a second pass necessary.
	print STDERR "Starting first BLAST pass for $Fld[0] - $Fld[1] at ";
	system("date");
	open FHR, "$blastall -C3 -F\"m S\" -i $Fld[0] -d $Fld[1] -p blastp -v $Fld[3] -b $Fld[3] -M $matrix -z 5000000 -m7 | ./$blastParser $score_cutoff |";

	%theHits = ();
	while (<FHR>) {
		$aLine = $_;
		chomp ($aLine);
		@words = split (/\s/, $aLine);

		if (exists ($theHits {$words [0]})) {
			$theHits {$words [0]} = $theHits {$words [0]}." ".$words [1];
		}
		else {
			$theHits {$words [0]} = $words [1];
		}

	}
	close (FHR);

	$tmpdir = ".";   # May be slightly (5%) faster using the RAM disk "/dev/shm".
	$tmpi = "$tmpdir/". $Fld[4] . "tmpi";
	$tmpd = "$tmpdir/". $Fld[4] . "tmpd";

	# Do second pass with compositional adjustment off to get full-length alignments.
	print STDERR "Starting second BLAST pass for $Fld[0] - $Fld[1] at ";
	system("date");
	unlink "$Fld[4]";
	$count = 0;
	$incr = 0;
	$size = keys %theHits;
	foreach $aQuery (sort keys % theHits) {
		# Create single-query file
		open (FHT, ">$tmpi");
		print FHT ">$aQuery\n".$sequencesA {">$aQuery"}."\n";
		close (FHT);
	  # Create mini-database of hit sequences
		open (FHT, ">$tmpd");
		foreach $aHit (split (/\s/, $theHits {$aQuery})) {
			print FHT ">$aHit\n".$sequencesB {">$aHit"}."\n";
		}
		close (FHT);
		# Run Blast and add to output
		$count++;
		if ($count/$size > $incr) {
		    print "Done ".100*$incr."% at ";
		    system("date");
		    $incr+=0.1;
		}
		system ("$formatdb -i $tmpd");
		system ("$blastall -C0 -FF -i $tmpi -d $tmpd -p blastp -v $Fld[3] -b $Fld[3] -M $matrix -z 5000000 -m7 | ./$blastParser $score_cutoff >> $Fld[4]");
	}
	unlink "$tmpi", "$tmpd", "formatdb.log", "$tmpd.phr", "$tmpd.pin", "$tmpd.psq";
}

sub do_diamond {
  use lib './';
  require "$diamondParser";
  my @Fld = @_;
  system("touch $Fld[4]");
  print STDERR "Formatting Diamond databases for $Fld[0] and $Fld[1]\n";
  system ("$diamond makedb --in $Fld[1] -d $Fld[1].db --quiet");
  if ($two_passes) {
    do_diamond_2pass(@_);
  } else {
    do_diamond_1pass(@_);
  }
}

sub do_diamond_1pass {
  my @Fld = @_;
  use lib './';
  require "$diamondParser";
  print STDERR "Starting Diamond searches...\n";
	if ($cores>0){
		system ("$diamond blastp -d $Fld[1].db.dmnd -q $Fld[0] -o diamond$Fld[0]$Fld[1] -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0 -p $cores");
	}else{
		system ("$diamond blastp -d $Fld[1].db.dmnd -q $Fld[0] -o diamond$Fld[0]$Fld[1] -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0");
	}
  diamondParser::parse("diamond$Fld[0]$Fld[1]", $score_cutoff, "$Fld[0]-$Fld[1]", 0,0);
  unlink "$Fld[1].db.dmnd", "diamond$Fld[0]$Fld[1]";
}

sub do_diamond_2pass {
  use lib './';
  require "$diamondParser";
  my @Fld = @_;
  %sequencesA = ();
  %sequencesB = ();

  open (FHA, $Fld [0]);
  while (<FHA>) {
		$aLine = $_;
		chomp ($aLine);
		$seq = "";
		if ($aLine =~ />/) {
			$header = $aLine;
			$sequencesA {$header} = "";
		}
		else {
			$sequencesA {$header} = $sequencesA {$header}.$aLine;
		}
  }
  close (FHA);

  open (FHB, $Fld [1]);
  while (<FHB>) {
		$aLine = $_;
		chomp ($aLine);
		$seq = "";
		if ($aLine =~ />/) {
			$header = $aLine;
			$sequencesB {$header} = "";
		}
		else {
			$sequencesB {$header} = $sequencesB {$header}.$aLine;
		}
  }
  close (FHB);
  my $tmpResults="tmp$Fld[0]-$Fld[1]";
  print STDERR "Starting first Diamond pass for $Fld[0] - $Fld[1] at ";
  system("date");
	if ($cores>0){
		system ("$diamond blastp -d $Fld[1].db.dmnd -q $Fld[0] -o tmp$Fld[0]$Fld[1] -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 4 --max-hsps 0 -p $cores");
	}else{
		system ("$diamond blastp -d $Fld[1].db.dmnd -q $Fld[0] -o tmp$Fld[0]$Fld[1] -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 4 --max-hsps 0");
	}
  diamondParser::parse("tmp$Fld[0]$Fld[1]", $score_cutoff, $tmpResults, 0,0);

  %theHits = ();
  open (my $tmpResults, "<", $tmpResults) or die "Could not open output file: $tmpResults \n ";
  while (<$tmpResults>) {
		$aLine = $_;
		chomp ($aLine);
		@words = split (/\t/, $aLine);
		if (exists ($theHits {$words [0]})) {
			$theHits {$words [0]} = $theHits {$words [0]}." ".$words [1];
		}
		else {
			$theHits {$words [0]} = $words [1];
		}
  }
  #print("\nNumber of hits: ", %theHits, "\n");
  $tmpdir = ".";   # May be slightly (5%) faster using the RAM disk "/dev/shm".
  $tmpi = $Fld[4] . "tmpi";
  $tmpd = $Fld[4] . "tmpd";

  # Do second pass with compositional adjustment off to get full-length alignments.
  print STDERR "Starting second Diamond pass for $Fld[0] - $Fld[1] at ";
  system("date");
  unlink "$Fld[4]";
  $count = 0;
  $incr = 0;
  $size = keys %theHits;
  foreach $aQuery (sort keys % theHits) {
		# Create single-query file
		open (FHT, ">$tmpi");
		print FHT ">$aQuery\n".$sequencesA {">$aQuery"}."\n";
		close (FHT);
	        # Create mini-database of hit sequences
		open (FHT, ">$tmpd");
		foreach $aHit (split (/\s/, $theHits {$aQuery})) {
			print FHT ">$aHit\n".$sequencesB {">$aHit"}."\n";
		}
		close (FHT);
		# Run Blast and add to output
		$count++;
		if ($count/$size > $incr) {
		    print "Done ".100*$incr."% at ";
		    system("date");
		    $incr+=0.1;
		}
  	system ("$diamond makedb --in $tmpd -d $tmpd --quiet");
		if ($cores>0){
			system ("$diamond blastp -d $tmpd.dmnd -q $tmpi -o $tmpd$$tmpi -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0 --threads $cores");
		}else{
			printf "Using all available cores";
			system ("$diamond blastp -d $tmpd.dmnd -q $tmpi -o $tmpd$$tmpi -f 6 qseqid sseqid bitscore qlen slen qstart qend sstart send --$sensitivity --quiet --comp-based-stats 0 --max-hsps 0");
		}
  	diamondParser::parse("$tmpd$$tmpi", $score_cutoff, $Fld[4], 1,0);
	  }
	  unlink "$Fld[1].db.dmnd", "$Fld[0].db.dmnd", "$Fld[1]-$Fld[0]tmpd.dmnd", "$Fld[0]-$Fld[1]tmpd.dmnd", "tmp$Fld[0]-$Fld[1]", "tmp$Fld[0]$Fld[1]", "tmp$Fld[1]$Fld[0]";
	  unlink "$tmpi", "$tmpd", "formatdb.log", "$tmpd.phr", "$tmpd.pin", "$tmpd.psq";
	}
}


#   Date                                 Modification
# --------          ---------------------------------------------------
#
# 2006-04-02 [1.36] - Changed score cutoff 50 to 0 for blast2faa.pl.
#                   Reason: after a cluster merger a score can be less than the cutoff (50)
#                   which will remove the sequence in blast2faa.pl.  The bootstrapping will
#                   then fail.
#                   - Fixed bug with index variable in bootstrap routine.
#
# 2006-06-01 [2.0]  - Fixed bug in blast_parser.pl: fields 7 and 8 were swapped,
#                   it was supposed to print match_area before HSP_length.
#                   - Fixed bug in blastall call: -v param was wrong for the A-B
#                   and B-A comparisons.
#                   -
#                   - Changed "cluster" to "group" consistently in output.
#                   - Changed "main ortholog" to "seed ortholog" in output.
#                   - Replace U -> X before running seqstat.jar, otherwise it crashes.
# 2006-08-04 [2.0]  - In bootstrap subroutine, replace U with X, otherwise seqstat
#                       will crash as this is not in the matrix (should fix this in seqstat)
# 2006-08-04 [2.1]  - Changed to writing default output to file.
#                   - Added options to run blast only.
#                   - Fixed some file closing bugs.
# 2007-12-14 [3.0]  - Sped up sorting routines (by Isabella).
#                   - New XML-based blast_parser.
#                   - New seqstat.jar to handle u and U.
#                   - Modified overlap criterion for rejecting matches.  Now it agrees with the paper.
# 2009-04-01 [4.0]  - Further modification of overlap criteria (require that they are met for both query and subject).
#		    - Changed bit score cutoff to 40, which is suitable for compositionally adjusted BLAST.
#		    - Added in 2-pass algorithm.
# 2009-06-11 [4.0]  - Moved blasting out to subroutine.
#		    - Changed blasting in bootstrap subroutine to use unconditional score matrix adjustment and SEG hard masking,
#		      to be the same as first step of 2-pass blast.
# 2009-06-17 [4.0]  - Compensated a Blast "bug" that sometimes gives a self-match lower score than a non-identical match.
#                      This can happen with score matrix adjustment and can lead to missed orthologs.
# 2009-08-18 [4.0]  - Consolidated Blast filtering parameters for 2-pass (-C3 -F\"m S\"; -C0 -FF)
# 2019-02-05 [4.2]  - Renamed output "Output" to "Stats"
#                   - Modifications for all-vs-all so that blast is not run if output exists.
#                   - Modifications for allowing parallel runs in the same directory, to use all CPUs.
#                   - Added sorted running of blast pass 2 for reproducibility.
#                   - Fixed bug in blast_parser.pl that produced an error if first query found no match.
# 2021-06-15 [5.0]  - Updated to work with DIAMOND for sequence alignment
