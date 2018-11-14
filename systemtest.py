import os, subprocess, platform

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
local_blast_path = os.path.join(cDir, 'dep-win64', 'ncbi-blast-2.2.26')

def test_java(output_msg):
    JAVA_OK = False
    try:
        out = subprocess.call(['java', '-version'])
        JAVA_OK = True
    except (OSError):
        output_msg.append('MetaDraft requires a working JRE, see README.md for details.')
    return JAVA_OK, output_msg

def test_perl(output_msg):
    PERL_OK = False
    try:
        out = subprocess.call(['perl', '-v'])
        PERL_OK = True
    except (OSError):
        output_msg.append('MetaDraft requires a working Perl installation, see README.md for details.')
    return PERL_OK, output_msg

def test_perl_xml(output_msg):
    PERL_XML_OK = False
    p_script = """\
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
    """
    try:
        PF = file('_test.pl', 'w')
        PF.write(p_script)
        PF.close()
        out = int(subprocess.call(['perl', '_test.pl']))
        if out:
            raise OSError
        PERL_XML_OK = True
    except (OSError):
        output_msg.append('MetaDraft requires Perl has the XML::Parser package installed, see README.md for details.')
    return PERL_XML_OK, output_msg

def test_blast(output_msg):
    BLAST_OK = BLAST_HAVE_LOCAL = False
    pth = ''
    try:
        out = subprocess.call(['formatdb'])
        #out = subprocess.call(['blastall'])
        BLAST_OK = True
    except (OSError):
        # test to see if we can use a built-in binary
        if os.name == 'nt' and platform.machine().endswith('64'):
            for pth in ['PATH', 'path', 'Path']:
                if pth in os.environ:
                    os.environ[pth] = os.environ[pth] + ';' + local_blast_path
                    break
            try:
                out = subprocess.call(['formatdb'])
                #out = subprocess.call(['blastall', '--help'])
                BLAST_OK = True
                BLAST_HAVE_LOCAL = True
            except (OSError):
                output_msg.append('MetaDraft requires a working NCBI BLAST2 in the path, see README.md for details.')
        else:
            output_msg.append('MetaDraft requires a working NCBI BLAST2 in the path, see README.md for details.')
    if os.path.exists(os.path.join(cDir, 'formatdb.log')):
        os.remove(os.path.join(cDir, 'formatdb.log'))
    if os.path.exists(os.path.join(cDir, '_test.pl')):
        os.remove(os.path.join(cDir, '_test.pl'))

    return BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg

def print_test_results(JAVA_OK, PERL_OK, PERL_XML_OK, BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg):
    print('\n\nSystem check results:\n=====================')
    print('Java test passed: {}'.format(JAVA_OK))
    print('Perl test passed: {}'.format(PERL_OK))
    print('Perl XML test passed: {}'.format(PERL_XML_OK))
    print('BLAST test passed: {}'.format(BLAST_OK))
    if not BLAST_OK:
        print('BLAST can use built-in test passed: {}'.format(BLAST_HAVE_LOCAL))

    if len(output_msg) > 0:
        print('\n\nSuggestions:\n============')
        for l in output_msg:
            print(l)
        print('\n')
    if BLAST_HAVE_LOCAL:
        print("MetaDraft requires NCBI BLAST and can make use of it's own distribution. Please consider adding \'{}\' to your local '{}' and see README.md for details.".format(local_blast_path, pth))

if __name__ == '__main__':
    output_msg = []
    JAVA_OK, output_msg = test_java(output_msg)
    PERL_OK, output_msg = test_perl(output_msg)
    PERL_XML_OK, output_msg = test_perl_xml(output_msg)
    BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg = test_blast(output_msg)
    print_test_results(JAVA_OK, PERL_OK, PERL_XML_OK, BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg)
    if JAVA_OK and PERL_OK and PERL_XML_OK and ( BLAST_OK or BLAST_HAVE_LOCAL ):
        os.sys.exit(0)
    else:
        os.sys.exit(1)


