"""
MetaToolkit: MetaDraft
======================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit: MetaDraft (https://github.com/SystemsBioinformatics/cbmpy-metadraft)
Copyright (C) 2016-2018 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands

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
    if RELEASE_STATUS == 2:
        widget_splash.showMessage("Ver {}-({}) beta\nAuthor: Brett G. Olivier PhD\n(c) Brett G. Olivier, Amsterdam, 2017-2018.\nSee Help - About for more details.".format(metadraft_version, cbmpy.__version__), alignment=Qt.AlignBottom)
        time.sleep(3)
    else:
        widget_splash.showMessage("Ver {}-({}) beta\n(c) Brett G. Olivier, Amsterdam, 2017-2018.\nSee Help - About for more details.".format(metadraft_version, cbmpy.__version__), alignment=Qt.AlignBottom)
        if RELEASE_STATUS == 1:
            time.sleep(2)
    ex = MetaDraftApp()
    widget_splash.finish(ex)
    sys.exit(app.exec_())


if __name__ == '__main__':
    import os
    import systemtest

    output_msg = []
    BLAST_OK, BLAST_HAVE_LOCAL, pth, output_msg = systemtest.test_blast(output_msg)
    if not ( BLAST_OK or BLAST_HAVE_LOCAL ):
        print(output_msg[0])
        os.sys.exit(1)
    elif BLAST_HAVE_LOCAL:
        print("\nMetaDraft requires NCBI BLAST but can make use of it's own distribution. I have set the PATH for you but please consider adding \'{}\' to your local '{}' environment variable to remove this message. Please see README.md for details.\n".format(systemtest.local_blast_path, pth))

    import libpython.qtmetadraft
    from libpython.qtmetadraft import *
    print(API_VERSION)
    __version__ = metadraft_version
    run_metadraft()


