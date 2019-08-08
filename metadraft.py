"""
MetaToolkit: MetaDraft
======================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit: MetaDraft (https://github.com/SystemsBioinformatics/cbmpy-metadraft)
Copyright (C) 2016-2019 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands

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
    app = QApplication(sys.argv)
    widget_splash = QSplashScreen(QPixmap("images/metatoolkit1-03.jpg"))
    widget_splash.show()
    widget_splash.showMessage("Ver {}-({})\nAuthor: Brett G. Olivier\n(c) Systems Bioinformatics, VU University Amsterdam, Amsterdam, 2017-2019.\nSee Help - About for more details.".format(metadraft_version, cbmpy.__version__), alignment=Qt.AlignBottom)
    time.sleep(2)
    ex = MetaDraftApp()
    widget_splash.finish(ex)
    sys.exit(app.exec_())


if __name__ == '__main__':
    import os, json, platform
    import systemtest

    F = open('_metadraft.cfg', 'r')
    config = json.load(F)
    F.close()
    if platform.architecture() == ('64bit', 'WindowsPE') and not config['system']['have_blas2']:
        output_msg = []
        BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg = systemtest.test_blast(output_msg)
        if not ( BLAST_OK or BLAST_HAVE_LOCAL ):
            print(output_msg[0])
            os.sys.exit(1)
        elif BLAST_OK and not BLAST_HAVE_LOCAL:
            config['system']['have_blas2'] = True
            F = open('_metadraft.cfg', 'w')
            json.dump(config, F)
            F.close()
        elif BLAST_HAVE_LOCAL:
            print("\nMetaDraft requires NCBI BLAST but can make use of it's own distribution. I have set the PATH for you but please consider adding \'{}\' to your local '{}' environment variable to remove this message. Please see README.md for details.\n".format(systemtest.local_blast_path, pth))

    import libpython.qtmetadraft
    from libpython.qtmetadraft import *
    print(API_VERSION)
    __version__ = metadraft_version
    run_metadraft()


