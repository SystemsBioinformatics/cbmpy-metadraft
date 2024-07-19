my $rc = 0;
    $rc = eval
    {
      require XML::Parser;
      XML::Parser->import();
      1;
    };
    if ($rc){
        exit 0
        } else {
        exit 1
        }