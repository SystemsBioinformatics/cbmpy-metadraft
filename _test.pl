    my $rc = 0;
    $rc = eval
    {
      require XML::Parser;
      XML::Parser->import();
      1;
    };
    if ($rc){
        print "happy";
        exit 0
        } else {
        print "sad";
        exit 1
        }
    