"""
MetaToolkit: MetaDraft
======================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit: MetaDraft (https://github.com/SystemsBioinformatics/cbmpy-metadraft)
Copyright (C) 2016-2024 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier

"""

import os, subprocess, platform

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
LOCAL_BLASTWIN_PATH = os.path.join(cDir, 'dep-win64', 'ncbi-blast-2.2.26')
LOCAL_DIAMOND_PATH =  os.path.join(cDir, 'diamond')
LOCAL_INPARANOID_PATH =  os.path.join(cDir, 'inparanoid')

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
        print "happy";
        exit 0
        } else {
        print "sad";
        exit 1
        }
    """
    try:
        PF = open('_test.pl', 'w')
        PF.write(p_script)
        PF.close()
        out = int(subprocess.call(['perl', '_test.pl']))
        print(out)
        if out:
            raise OSError
        PERL_XML_OK = True
    except (OSError):
        output_msg.append('MetaDraft requires Perl has the XML::Parser package installed, see README.md for details.')

    # if os.path.exists(os.path.join(cDir, '_test.pl')):
        # os.remove(os.path.join(cDir, '_test.pl'))

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
                    os.environ[pth] = os.environ[pth] + ';' + LOCAL_BLASTWIN_PATH
                    break
            try:
                out = subprocess.call(['formatdb'])
                #out = subprocess.call(['blastall', '--help'])
                BLAST_OK = True
                BLAST_HAVE_LOCAL = True
            except (OSError):
                output_msg.append('Ignore: MetaDraft requires a working NCBI BLAST2 in the path, see README.md for details.')
        else:
            output_msg.append('Ignore: MetaDraft requires a working NCBI BLAST2 in the path, see README.md for details.')
    if os.path.exists(os.path.join(cDir, 'formatdb.log')):
        os.remove(os.path.join(cDir, 'formatdb.log'))

    output_msg.append('BLAST support has been deprecated and is no longer necessary. MetaDraft now uses InParanoid-DIAMOND, see README.md for details. The output of this test can be ignored.')

    return BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg

def test_python_dependencies(output_msg):
    PYTHON_DEP_OK = True
    HAVE_QT4 = HAVE_QT5 = HAVE_QT6 = False
    try:
        import PyQt4
        HAVE_QT4 = True
    except ImportError:
        pass
        #output_msg.append('PyQt4 not found, see README.md for details.')
    try:
        import PyQt5
        HAVE_QT5 = True
    except ImportError:
        pass
        #output_msg.append('PyQt5 not found, see README.md for details.')
    try:
        import PyQt6
        HAVE_QT6 = True
    except ImportError:
        pass
        #output_msg.append('PyQt5 not found, see README.md for details.')

    if HAVE_QT6:
        output_msg.append('PyQt6 found, excellent, this is the recommended version.')
    if HAVE_QT5:
        output_msg.append('PyQt5 found, MetaDraft supports this but will move to PyQT6 in the future.')
    if HAVE_QT4 and not HAVE_QT5 and not HAVE_QT6:
        output_msg.append('PyQt4 support has been deprecated, please install PyQt5 or PyQt6.')
        PYTHON_DEP_OK = False
    if not HAVE_QT5 and not HAVE_QT6:
        output_msg.append('Supported PyQt not found, please install PyQt5 or PyQt6 (recommended) before running MetaDraft, see README.md for details.')
        PYTHON_DEP_OK = False

    try:
        import libsbml
        import cbmpy
        import Bio
        import xlrd, xlwt, xlsxwriter
    except ImportError:
        output_msg.append('MetaDraft requires CBMPy, libSBML, xlsxwriter and BioPython to be installed. Please install before running MetaDraft, see README.md for details.')
        PYTHON_DEP_OK = False

    return PYTHON_DEP_OK, output_msg


def test_diamond(output_msg):
    DIAMOND_OK = DIAMOND_HAVE_LOCAL = False
    pth = ''
    try:
        out = subprocess.call(['diamond'])
        DIAMOND_OK = True
    except (OSError):
        # test to see if we can use a built-in binary
        if os.name == 'posix':
            for pth in ['PATH', 'path', 'Path']:
                if pth in os.environ:
                    os.environ[pth] = os.environ[pth] + ':' + LOCAL_DIAMOND_PATH
                    print(os.environ[pth])
                    break
            try:
                out = subprocess.call(['diamond'])
                DIAMOND_OK = True
                DIAMOND_HAVE_LOCAL = True
            except (OSError):
                output_msg.append('MetaDraft requires a working "diamond" in the path, see README.md for details.')
        elif os.name == 'nt':
            for pth in ['PATH', 'path', 'Path']:
                if pth in os.environ:
                    os.environ[pth] = os.environ[pth] + ';' + LOCAL_DIAMOND_PATH
                    break
            try:
                out = subprocess.call(['diamond.exe'])
                DIAMOND_OK = True
                DIAMOND_HAVE_LOCAL = True
            except (OSError):
                output_msg.append('MetaDraft requires a working "diamond" in the path, see README.md for details.')
        else:
            output_msg.append('MetaDraft requires a working "diamond" in the path, see README.md for details.')

    return DIAMOND_OK, DIAMOND_HAVE_LOCAL, pth, output_msg


def print_test_results(*args):
    print('\n\nSystem check results:\n=====================')
    print('Java test passed: {}'.format(JAVA_OK))
    print('Perl test passed: {}'.format(PERL_OK))
    print('Perl XML test passed: {}'.format(PERL_XML_OK))
    print('Python dependency test passed: {}'.format(PYTHON_DEP_OK))
    print('BLAST test passed (deprecated): {}'.format(BLAST_OK))
    print('DIAMOND test passed: {}'.format(DIAMOND_OK))
    if DIAMOND_OK:
        print('DIAMOND using local: {}'.format(DIAMOND_HAVE_LOCAL))

    if len(output_msg) > 0:
        print('\n\nSuggestions:\n============')
        for l in output_msg:
            print(l)
        print('\n')
    if BLAST_HAVE_LOCAL:
        print("MetaDraft requires NCBI BLAST and can make use of it's own distribution. Please consider adding \'{}\' to your local '{}' and see README.md for details.".format(LOCAL_BLASTWIN_PATH, pth))

    if os.name == 'nt' and not JAVA_OK:
        print('Make sure OpenJDk is in your path, for example: path=%path%;c:\DevTools\jdk-23.0.2\bin')
    elif os.name == 'nt' and not PERL_OK:
        print('Make sure perl is in your path, for example: path=%path%;c:\DevTools\strawberry\c\bin;c:\DevTools\strawberry\perl\bin')



if __name__ == '__main__':
    output_msg = []
    JAVA_OK, output_msg = test_java(output_msg)
    PERL_OK, output_msg = test_perl(output_msg)
    PERL_XML_OK, output_msg = test_perl_xml(output_msg)
    PYTHON_DEP_OK, output_msg = test_python_dependencies(output_msg)
    BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg = test_blast(output_msg)
    DIAMOND_OK, DIAMOND_HAVE_LOCAL, pth, output_msg = test_diamond(output_msg)


    print_test_results(JAVA_OK, PERL_OK, PERL_XML_OK, BLAST_OK, BLAST_HAVE_LOCAL, PYTHON_DEP_OK, DIAMOND_OK, DIAMOND_HAVE_LOCAL, \
                       pth, output_msg)
    if JAVA_OK and PERL_OK and PERL_XML_OK and ( DIAMOND_OK or DIAMOND_HAVE_LOCAL ) and PYTHON_DEP_OK:
        print('\nCongratulations you are ready to run MetaDraft! (python metadraft.py)\n')
        os.sys.exit(0)
    else:
        print('\nYou need to install a few more things before you are ready to run MetaDraft.\n')
        os.sys.exit(1)


