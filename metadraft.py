"""
MetaToolkit: MetaDraft
======================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit: MetaDraft (https://github.com/SystemsBioinformatics/cbmpy-metadraft)
Copyright (C) 2016-2023 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands

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

Author: Brett G. Olivier PhD
Contact email: b.g.olivier@vu.nl

"""
from __future__ import division, print_function
from __future__ import absolute_import


def run_metadraft():
    """
    Launches the MetaDraft application with a splash screen.

    This function initializes a QApplication, displays a splash screen with an image and version information, and then runs the MetaDraft application. The version information includes MetaDraft version and CBMPy version. The splash screen remains visible for 2 seconds, after which the MetaDraft main application window is displayed.

    Args:
    None

    Returns:
    None

    Raises:
    sys.exit: Exits the application normally after the QApplication event loop ends.
    """
    app = QApplication(sys.argv)
    widget_splash = QSplashScreen(QPixmap("images/binaries.jpg"))
    widget_splash.show()
    
    widget_splash.showMessage(
        "Ver {}-({})\nAuthor: Brett G. Olivier\n(c) A-LIFE, VU University Amsterdam, Amsterdam, 2017-2025.\nSee Help - About for more details.".format(
            metadraft_version, cbmpy.__version__
        ),
        alignment=QtCore.Qt.AlignmentFlag.AlignBottom,
    )
    time.sleep(0.5)
    ex = MetaDraftApp()
    widget_splash.finish(ex)
    sys.exit(app.exec())


if __name__ == '__main__':
    import os, json, platform
    import systemtest

    F = open('_metadraft.cfg', 'r')
    config = json.load(F)
    F.close()
    #if (
        #platform.architecture() == ('64bit', 'WindowsPE')
        #and not config['system']['have_blas2']
    #):
        #output_msg = []
        #BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg = systemtest.test_blast(output_msg)
        #if not (BLAST_OK or BLAST_HAVE_LOCAL):
            #print(output_msg[0])
            #os.sys.exit(1)
        #elif BLAST_OK and not BLAST_HAVE_LOCAL:
            #config['system']['have_blas2'] = True
            #F = open('_metadraft.cfg', 'w')
            #json.dump(config, F)
            #F.close()
        #elif BLAST_HAVE_LOCAL:
            #print(
                #"\nMetaDraft requires NCBI BLAST but can make use of it's own distribution. I have set the PATH for you but please consider adding \'{}\' to your local '{}' environment variable to remove this message. Please see README.md for details.\n".format(
                    #systemtest.LOCAL_BLASTWIN_PATH, pth
                #)
            #)


    DIAMOND_OK, DIAMOND_HAVE_LOCAL, pth, output_msg = systemtest.test_diamond([])
    print('DIAMOND_OK',  DIAMOND_OK)
    print('DIAMOND_HAVE_LOCAL', DIAMOND_HAVE_LOCAL)
    print('pth', pth)
    if not DIAMOND_OK:
        print('InParanoid-DIAMOND is needed for this version of metadraft')
        print('DIAMOND_HAVE_LOCAL', DIAMOND_HAVE_LOCAL)
        os.sys.exit(2)

    # import libpython.qtmetadraft
    from libpython.qtmetadraft import *

    __version__ = metadraft_version
    
    run_metadraft()
