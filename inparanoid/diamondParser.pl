package diamondParser;
use warnings;
use strict;

sub parse {
  my @input = @_;
  my $file1 = $input[0];
  my $score_cutoff = $input[1];
  my $output = $input[2];
  my $secondPass = $input[3];
  my $alignment_mode = $input[4];

  #remove the output file if it exists so it can be overwritten without risking appending to an older file
  if ($secondPass==0){
    if (-e $output) {
      unlink $output;
    }
  }


  # Open input and output files
  open (my $alignments, "<", $file1) or die "Can't find $file1 \n";
  open (my $out, ">>$output") or die "Could not open output file: $output \n ";

  # Read all lines
  my %hits=();
  my %hitsSequences=();
  while (my $line = <$alignments>) {
    chomp($line);
    my @words=split("\t", $line);

    my $queryName = $words[0];
    my $subjectName = $words[1];
    my $score = $words[2];
    my $queryLength = $words[3];
    my $subjectLength = $words[4];
    my $queryMatchingRegion = $words[6]-$words[5]+1;
    my $subjectMatchingRegion = $words[8]-$words[7]+1;
    my $querySubjectRegions = join("", "q:", $words[5], "-", $words[6], " h:", $words[7], "-", $words[8]);
    if ($alignment_mode==1) {
    		my $querySequence = $words[9];
    		my $subjectSequence = $words[10];
        my $btop = $words[11];
        # Getting alignment from btop
        my $alignedQ="";
        my $alignedS="";
        my $btopIndex=0;
        my $qIndex=0;
        my $sIndex=0;
        my $lastWasQ=0;
        my $longerThanOneCharNumber="";
        foreach my $char (split('', $btop)) {
          if( $char =~ /^\d*$/){
            if (length($btop)==$btopIndex+1 or substr($btop, $btopIndex+1, 1) =~ /\D/){
              $longerThanOneCharNumber=$longerThanOneCharNumber.$char;
              $alignedQ=$alignedQ.substr($querySequence,$qIndex,$longerThanOneCharNumber);
              $qIndex=$qIndex+$longerThanOneCharNumber;
              $alignedS=$alignedS.substr($subjectSequence,$sIndex,$longerThanOneCharNumber);
              $sIndex=$sIndex+$longerThanOneCharNumber;
              $longerThanOneCharNumber=""
            }else{
              $longerThanOneCharNumber=$longerThanOneCharNumber.$char;
            }
          }else{
            if ($lastWasQ==0){
              $alignedQ=$alignedQ.$char;
              $lastWasQ=1;
              if ($char ne "-"){
                $qIndex+=1;
              }
            }else{
              $alignedS=$alignedS.$char;
              $lastWasQ=0;
              if ($char ne "-"){
                $sIndex+=1;
              }
            }
          }
          $btopIndex+=1;
        }
    		if (exists $hitsSequences{"$queryName:$subjectName"}){
    			$hitsSequences{"$queryName:$subjectName"}[0]=$hitsSequences{"$queryName:$subjectName"}[0].$alignedQ;
    			$hitsSequences{"$queryName:$subjectName"}[1]=$hitsSequences{"$queryName:$subjectName"}[1].$alignedS;
    		}else{
    			$hitsSequences{"$queryName:$subjectName"}=[$alignedQ,$alignedS];
    		}
	   }

    # Do not include hits with lower score than the cutoff
    my $overlapping = 0;
    # If the hit already exist, they need to be merged
    if (exists $hits{"$queryName:$subjectName"}){
        my @existingHit=@{$hits{"$queryName:$subjectName"}};
        my @existingAlign=@existingHit[9..$#existingHit];
        my $rlowestQ = $words[5];
        my $rhighestQ = $words[6];
        my $rlowestS = $words[7];
        my $rhighestS = $words[8];
	my $rlengthQ=$rhighestQ-$rlowestQ+1;
	my $rlengthS=$rhighestS-$rlowestS+1;

	my $lowestQ = $rlowestQ;
        my $highetsQ = $rhighestQ;
        my $lowestS = $rlowestS;
        my $highetsS = $rhighestS;

        foreach my $align (@existingAlign){

	my $otherlowestQ = (split /-/, ((split /:/, $align)[1]))[0];
        my $otherhighestQ =(split / /, (split /-/, ((split /:/, $align)[1]))[1])[0];
        my $otherlowestS =(split /-/, ((split /:/, $align)[2]))[0];
        my $otherhighestS = (split /-/, ((split /:/, $align)[2]))[1];
	my $otherlengthQ=$otherhighestQ-$otherlowestQ+1;
	my $otherlengthS=$otherhighestS-$otherlowestS+1;

	#Hits are allowed to overlap with 5% of the shortest sequence
 	my $shortest_lengthQ = ($otherlengthQ < $rlengthQ)?$otherlengthQ:$rlengthQ;
	my $shortest_lengthS = ($otherlengthS < $rlengthS)?$otherlengthS:$rlengthS;

	my $fivePercentOfShortestQ = $shortest_lengthQ*0.05;
        my $fivePercentOfShortestS = $shortest_lengthS*0.05;

          if($lowestQ>(split /-/, ((split /:/, $align)[1]))[0]){
            $lowestQ=(split /-/, ((split /:/, $align)[1]))[0];
          }
          if($highetsQ<(split / /, (split /-/, ((split /:/, $align)[1]))[1])[0]){
            $highetsQ=(split / /, (split /-/, ((split /:/, $align)[1]))[1])[0];
          }
          if($lowestS>(split /-/, ((split /:/, $align)[2]))[0]){
            $lowestS=(split /-/, ((split /:/, $align)[2]))[0];
          }
          if($highetsS<(split /-/, ((split /:/, $align)[2]))[1]){
            $highetsS=(split /-/, ((split /:/, $align)[2]))[1];
          }

	 if (($rlowestQ <= $otherhighestQ-$fivePercentOfShortestQ && $rhighestQ >= $otherhighestQ+$fivePercentOfShortestQ) || ($rhighestQ >= $otherlowestQ+$fivePercentOfShortestQ && $rlowestQ <= $otherlowestQ-$fivePercentOfShortestQ) || ($rlowestQ>= $otherlowestQ-$fivePercentOfShortestQ && $rhighestQ <= $otherhighestQ+$fivePercentOfShortestQ)){
                    $overlapping=1;
}
                if (($rlowestS <= $otherhighestS-$fivePercentOfShortestS && $rhighestS >= $otherhighestS+$fivePercentOfShortestS) || ($rhighestS >= $otherlowestS+$fivePercentOfShortestS && $rlowestS <= $otherlowestS-$fivePercentOfShortestS) || ($rlowestS >= $otherlowestS-$fivePercentOfShortestS && $rhighestS <= $otherhighestS+$fivePercentOfShortestS)){
                   $overlapping=1;
}
                if (($rlowestQ > $otherhighestQ-$fivePercentOfShortestQ && $rlowestS < $otherhighestS-$fivePercentOfShortestS) || ($rhighestQ < $otherlowestQ+$fivePercentOfShortestQ && $rhighestS > $otherlowestS+$fivePercentOfShortestS)){
                    $overlapping=1;
}
      }
	     if ($overlapping==0){
        $hits{"$queryName:$subjectName"}=[$queryName, $subjectName, $score+$existingHit[2], $queryLength, $subjectLength, $highetsQ-$lowestQ+1, $highetsS-$lowestS+1, $existingHit[7]+$queryMatchingRegion, $existingHit[8]+$subjectMatchingRegion, @existingHit[9..$#existingHit], $querySubjectRegions];
      }
      }else{
        $hits{"$queryName:$subjectName"}=[$queryName, $subjectName, $score, $queryLength, $subjectLength, $queryMatchingRegion, $subjectMatchingRegion, $queryMatchingRegion, $subjectMatchingRegion, $querySubjectRegions];
      }
  }
  # Sort by name and score, and print to file
  my %hitScoreHash=();
  foreach my $h ( keys %hits){
    my $sortName=join(":", (split /:/, $h)[0], @{$hits{$h}}[2], (split /:/, $h)[1]);
    $hitScoreHash{$sortName}=join("\t", @{$hits{$h}});
  }
  foreach my $h (sort char_then_num keys %hitScoreHash){
    if ($alignment_mode==1) {
      print $out ">";
    }
    print $out $hitScoreHash{$h},"\n";
    if ($alignment_mode==1) {
      my $hitNames = join(":", (split /:/, $h)[0], (split /:/, $h)[2]);
    	print $out $hitsSequences{$hitNames}[0],"\n";
    	print $out $hitsSequences{$hitNames}[1],"\n";
    }
  }
}


# Function to sort by character and then score
sub char_then_num {
    my ($a_char, $a_num) = split ':', $a;
    my ($b_char, $b_num) = split ':', $b;
    return $a_char cmp $b_char
             ||
        $b_num <=> $a_num;
}


# due to how perl works the last line in a module must return true so don't remove the random 1 here
1;
