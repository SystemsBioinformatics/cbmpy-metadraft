"""
MetaToolkit: MetaDraft module
=============================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit:MetaDraft (http://cbmpy.sourceforge.net)
Copyright (C) 2015-2017 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Contact email: bgoli@users.sourceforge.net

"""

## http://qt-project.org/forums/viewthread/26126
## http://stackoverflow.com/questions/2304199/how-to-sort-a-qtablewidget-with-my-own-code
## http://qt-project.org/doc/qt-4.8/qtablewidgetitem.html

import os, sys, random, json, itertools, datetime, re, logging, webbrowser, copy
import zipfile, json, shutil, subprocess, math, time, threading, pprint, stat, csv

try:
    type(reduce)
except NameError:
    from functools import reduce

import sip
API_NAMES = ["QDate", "QDateTime", "QString", "QTextStream", "QTime", "QUrl", "QVariant"]
API_VERSION = 2
for name in API_NAMES:
    sip.setapi(name, API_VERSION)

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import numpy
import biotools
import biotoolsINP as bionoid
import metadraftreports
import cbmpy

try:
    import docx
    HAVE_DOCX = True
except ImportError:
    HAVE_DOCX = False


__version__ = '0.8.1'

HAVE_QT4 = False
HAVE_QT5 = False
try:
    import PyQt4
    from PyQt4 import (QtCore, QtGui, QtSvg)
    from PyQt4.QtCore import (Qt, QCoreApplication)
    from PyQt4.QtGui import (QMainWindow, QApplication, QPushButton, QWidget, QAction, QSplashScreen, QMdiArea)
    from PyQt4.QtGui import (QTabWidget, QVBoxLayout, QMdiSubWindow, QTableWidgetItem, QTextBrowser, QDockWidget)
    from PyQt4.QtGui import (QLabel, QTableWidget, QLineEdit, QComboBox, QGridLayout, QHBoxLayout, QButtonGroup)
    from PyQt4.QtGui import (QRadioButton, QFileDialog, QAbstractItemView, QMessageBox, QMenu, QSizePolicy)
    from PyQt4.QtGui import (QTextEdit, QListWidget, qApp, QStatusBar, QInputDialog)
    from PyQt4.QtGui import (QTreeView, QTreeWidget, QTreeWidgetItem, QProgressDialog)
    from PyQt4.QtGui import (QCursor, QBrush, QColor, QPalette, QPixmap, QFont)
    from PyQt4.QtCore import (pyqtSignal, pyqtSlot)
    HAVE_QT4 = True
    print('Using Qt4 - MetaMod is now also Qt5 compatible')
except ImportError as ex:
    pass

try:
    import PyQt5
    from PyQt5 import (QtCore, QtGui, QtSvg, QtWidgets)
    from PyQt5.QtCore import (Qt, QCoreApplication)
    from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, QAction, QSplashScreen, QMdiArea)
    from PyQt5.QtWidgets import (QTabWidget, QVBoxLayout, QMdiSubWindow, QTableWidgetItem, QTextBrowser, QDockWidget)
    from PyQt5.QtWidgets import (QLabel, QTableWidget, QLineEdit, QComboBox, QGridLayout, QHBoxLayout, QButtonGroup)
    from PyQt5.QtWidgets import (QRadioButton, QFileDialog, QAbstractItemView, QMessageBox, QMenu, QSizePolicy)
    from PyQt5.QtWidgets import (QTextEdit, QListWidget, qApp, QStatusBar, QInputDialog)
    from PyQt5.QtWidgets import (QTreeView, QTreeWidget, QTreeWidgetItem, QProgressDialog)
    from PyQt5.QtGui import (QCursor, QBrush, QColor, QPalette, QPixmap, QFont)
    from PyQt5.QtCore import (pyqtSignal, pyqtSlot)
    HAVE_QT5 = True
    print('Using Qt5 - the next generation')
except ImportError as ex:
    pass

if not HAVE_QT4 and not HAVE_QT5:
    print('Neither Qt4 nor Qt5 detected, please make sure PyQt is installed.')
    os.sys.exit(-1)

#import cbmpy as cbm

class NumberTableWidgetItem(QTableWidgetItem):
    def __lt__(self, other):
        if ( isinstance(other, QTableWidgetItem) ):
            my_value = self.data(Qt.EditRole)
            other_value = other.data(Qt.EditRole)
            try:
                my_value = float(str(my_value))
                my_ok = True
            except ValueError:
                my_ok = False
            try:
                other_value = float(str(other_value))
                other_ok = True
            except ValueError:
                other_ok = False

            if ( my_ok and other_ok ):
                return my_value < other_value

        return super(NumberTableWidgetItem, self).__lt__(other)

class NumberTableListLengthItem(QTableWidgetItem):
    def __lt__(self, other):
        if ( isinstance(other, QTableWidgetItem) ):
            my_value = len(self.data(Qt.EditRole).split(','))
            other_value = len(other.data(Qt.EditRole).split(','))
            #print(my_value, other_value)
            if ( my_value and my_value ):
                return my_value < other_value
        return super(NumberTableListLengthItem, self).__lt__(other)

class MyPopup(QWidget):
    def __init__(self):
        QWidget.__init__(self)

    def paintEvent(self, e):
        dc = QPainter(self)
        dc.drawLine(0, 0, 100, 100)
        dc.drawLine(100, 0, 0, 100)

class FileTreeView(QTreeView):
    def __init__(self, parent=None, rootpath=None):
        QTreeView.__init__(self)
        model = QFileSystemModel()
        model.setRootPath(rootpath)
        self.setModel(model)
        self.setRootIndex(model.index(rootpath))

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    http://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/

    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())


## 0: developer, 1: partner, 2: public
RELEASE_STATUS = 1

if RELEASE_STATUS > 0:
    __DEBUG__ = False
else:
    __DEBUG__ = True

class MetaDraftGUI(QWidget):
    appwindow = None
    stderr_logger = None
    stdout_logger= None
    log_files = os.path.join(cDir,'log')
    log_syslog = None
    cDir = cDir
    _tmpDir_ = os.path.join(cDir, 'tmp')
    __SYSLOG_ENABLED__ = True
    __DEBUG__ = __DEBUG__
    __DEL_BLAST_TMP__ = True
    _next_id_ = 1
    _DAT_MODELS = None
    _DAT_SEARCH_RES = None
    _DAT_G2REACT = None
    _DAT_LINK_DICT_ = None
    _CONFIG_ = None
    _SRC_TRG_MAP_ = None
    _DAT_NOTESDB_KEY_ = None

    NO_EXPORT_SEQ = True

    grid = None
    status_bar = None
    _colourCycler_ = None
    _colours_ = ((255, 255, 204), (204, 255, 204), (204, 255, 255), (204, 204, 255), (255, 204, 204),\
                 (204, 255, 153), (153, 255, 204), (102, 204, 255), (204, 153, 255), (255, 153, 204))
    _loading_ = None
    _active_tab_ = 0
    _last_tab_ = 0
    _active_tab_right_ = 0
    _last_tab_right_ = 0
    _tabpanel_idx_ = None
    _tabpanel_right_idx_ = None

    _wnotes_ = None
    notes = ''
    _wnotes_current_obj_ = None
    _wnotes_previous_obj_ = None

    _gene_selected_ids_ = None
    _gene_selected_map_ = None
    _reaction_selected_map_ = None
    _reaction_selected_ids_ = None
    _metab_selected_map_ = None
    selected_metabolites = None
    selected_reactions = None
    _gene_status_changed_ = False
    _metab_status_changed_ = False
    _session_actions_ = None
    _msg_box_ = None
    _tmp_cdir_template_ = None
    _widget_result_tree_rightclick_data_ = None
    _gene_score_limits_ = None
    _history_save_dir_ = cDir
    _history_open_dir_ = cDir

    model_name = 'sbml_model'
    default_compartment = 'cell'
    sbml_model = None
    sbml_save_dir = cDir
    sbml_file_name = ''
    excel_save_dir = cDir
    excel_file_name = ''
    data_input_files = cDir
    seqplus_files = os.path.join(cDir, 'lib_model')
    result_files = os.path.join(cDir, 'blast_results')
    blast_work_dir = os.path.join(cDir, 'data_blast')
    blast_tools = os.path.join(cDir, 'bin_base')
    metaproteome_files = os.path.join(cDir, 'lib_metaproteome')
    link_file = ''
    metaproteome_file = ''
    result_file = ''
    optimized_metaproteome = True
    _genedb_sqlcols_ = None
    _genedb_ = None
    _dummy_gene_obj_ = None
    _dbx_dir_ = os.path.join(cDir, 'dbx')
    _genedb_path_ = os.path.join(_dbx_dir_, '_metadraft_genedb.sql')
    _notesdb_path_ = os.path.join(_dbx_dir_, '_metadraft_notesdb.sql')
    config_file = '_metadraft.cfg'
    regex = None
    id_sep = '@'
    #gene_prefix = 'in_'
    gene_prefix = 'in_'
    reaction_table_loading = False
    gene_table_loading = False
    CURRENT_SELECTION_STATE = None
    REACT_LAST_CHECKED_ROW = 0
    GENE_LAST_CHECKED_ROW = 0
    METAB_LAST_CHECKED_ROW = 0

    REACT_LAST_SELECTED_ROW = 0
    GENE_LAST_SELECTED_ROW = 0
    METAB_LAST_SELECTED_ROW = 0

    GENE_SELECTION_STATE_CHANGE = False
    REACTION_SELECTION_STATE_CHANGE = False
    metadraft_rpt_footer = '<p><small>Generated by MetaDraft {} &copy; <a href="https://systemsbioinformatics.github.io/metatoolkit/">MetaToolkit</a>.</small></p>'.format(__version__)

    def __init__(self, appwindow):
        self.appwindow = appwindow
        self._loading_ = True

        super(MetaDraftGUI, self).__init__()
        #self.setAttribute(Qt.WA_DeleteOnClose)

        self._readConfig()
        ## debug
        if not self.__DEBUG__ and self._CONFIG_['system']['syslog']:
            self.__SYSLOG_ENABLED__ = True
            self.initSysLog()

        # create subdirectories for various things (they should exist)
        for f_ in [self.seqplus_files, self.result_files, self.blast_work_dir, self.blast_tools, self.metaproteome_files, self._tmpDir_]:
            if not os.path.exists(f_):
                os.makedirs(f_)

        # create default user result directory
        if not os.path.exists(os.path.join(self.result_files, self.func_getCurrentUser())):
            os.makedirs(os.path.join(self.result_files, self.func_getCurrentUser()))

        self.regex = {'GOterm' : re.compile('GO:\\d{7}')}

        # initialise the gene DB
        self._dummy_gene_obj_ = cbmpy.CBModel.Gene('dummy_gene')
        self._dummy_gene_obj_.setName('dummy_gene')
        self._dummy_gene_obj_.__setObjRef__(self._dummy_gene_obj_)
        self._dummy_gene_obj_.__global_id__ = {}
        def f(a, b, c):
            pass
        self._dummy_gene_obj_.__changeGlobalId__ = f
        # sync with metamod
        self._genedb_sqlcols_ = ['id TEXT PRIMARY KEY', 'pid TEXT', 'type TEXT', 'annotation TEXT', 'db_xref TEXT', 'notes TEXT', 'rdf TEXT', 'sbo TEXT', 'seq TEXT']
        self._notesdb_sqlcols_ = ['unixtime REAL PRIMARY KEY', 'user TEXT', 'model TEXT', 'id TEXT', 'notes TEXT', 'other TEXT']

        self._gene_score_limits_ = [1.0, 1.0]
        self._gene_selected_ids_ = []
        self._reaction_selected_map_ = {}
        self.selected_reactions = {}
        self._reaction_selected_ids_ = []
        self.selected_metabolites = []
        self._metab_selected_map_ = {}
        self.table_gene_cols = ['source', 'match', 'score', ' ', 'org']
        self.table_react_cols = ['reaction', 'name', ' ', 'org', 'genes', 'src']
        self.table_metab_cols = ['metabolite', 'name', 'org', 'fixed']
        self.grid = QGridLayout()
        self.grid.setSpacing(10)

        # create stuff
        self.createMenus()
        self.widgetTabPanel()
        self.widgetTabPanelRight()
        self.widgetBuildPanel()
        self.widgetDisplayReact()
        self.widgetBusy()

        # not needed with new applayout
        #spacer = QLabel(' ')
        #spacer.setMaximumHeight(20)
        #self.grid.addWidget(spacer, 0, 0, 1, 4)
        #self.grid.addWidget(spacer, 0, 4, 1, 3)

        # add tabs in order
        self._tabpanel_idx_ = {0 : 'Build',
                               1 : 'Genes',
                               2 : 'Reactions',
                               3 : 'Metabolites',
                               4 : 'Objective',
                               5 : 'Other'
                               }

        self.widgetTabPanel_add(self.build_panel, 'Build', setcurrent=True)
        self.grid.addWidget(self.widget_tabpanel,  1-1, 0, 4, 4)

        # add tabs in order
        self._tabpanel_right_idx_ = {0 : 'Information',
                                     1 : 'BuildReaction'
                                     }
        tab0R = self.widgetTabPanelRight_createtab(self.reactDisplay)
        self.widgetTabPanelRight_add(tab0R, 'Information', setcurrent=True)
        #tab1R = self.widgetTabPanelRight_createtab(QWidget())
        #self.widgetTabPanelRight_add(tab1R, 'BuildReaction', setcurrent=False)

        self.widget_tabpanel_right.setTabEnabled(1, False)
        #self.grid.addWidget(self.widget_tabpanel_right,  1-1, 4, 4, 3)
        self.grid.addWidget(self.widget_tabpanel_right,  1-1, 4, 3, 3)

        # Notes panel
        self._wnotes_ = QTextEdit()
        self._wnotes_.setMaximumHeight(60)
        self._wnotes_.setDisabled(True)
        self.grid.addWidget(self._wnotes_, 4-1, 4, 1, 3)

        self.setLayout(self.grid)

        self.status_bar = QStatusBar()
        self.status_bar.showMessage("Welcome to MetaDraft {}-{}".format(__version__, cbmpy.__version__))
        self.appwindow.setStatusBar(self.status_bar)

        # todo
        #self.widgetButtonPanel()

        # setup
        if RELEASE_STATUS == 2:
            self.appwindow.setWindowTitle('MetaDraft {} - created and developed by Brett G. Olivier PhD (b.g.olivier@vu.nl)'.format(__version__))
        else:
            self.appwindow.setWindowTitle('MetaDraft {} - Brett G. Olivier (b.g.olivier@vu.nl)'.format(__version__))
        self.initDBs()
        self._loading_ = False

    def _readConfig(self):
        if os.path.exists(os.path.join(self.cDir, self.config_file)):
            F = open(self.config_file, 'r')
            self._CONFIG_ = json.load(F)
            F.close()
        else:
            self._CONFIG_ = {}

    def _writeConfig(self):
        F = open(self.config_file, 'w')
        json.dump(self._CONFIG_, F)
        F.close()

    def initSysLog(self):
        if not os.path.exists(self.log_files):
            os.makedirs(self.log_files)
        self.log_syslog = os.path.join(self.log_files, 'syslog-{}'.format(time.strftime("%y-%m-%d")))
        logging.basicConfig(
            level = logging.DEBUG,
            format = '%(levelname)s:%(name)s:%(message)s',
            filename = self.log_syslog,
            filemode = 'a'
        )

        self.stdout_logger = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
        sys.stdout = self.stdout_logger

        self.stderr_logger = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
        sys.stderr = self.stderr_logger

    def initDBs(self):
        # GENEDB
        self._genedb_ = cbmpy.CBNetDB.DBTools()
        if not os.path.exists(os.path.join(self.cDir, self._genedb_path_)):
            print('INITDB: initialising geneDB')
            self._genedb_.connectSQLiteDB(self._genedb_path_)
            self._genedb_.createDBTable('GENES',  self._genedb_sqlcols_)
        else:
            print('INITDB: connecting geneDB')
            self._genedb_.connectSQLiteDB(self._genedb_path_)
        # NOTEDB
        self._notesdb_ = cbmpy.CBNetDB.DBTools()
        if not os.path.exists(os.path.join(self.cDir, self._notesdb_path_)):
            print('INITDB: initialising notesDB')
            self._notesdb_.connectSQLiteDB(self._notesdb_path_)
            self._notesdb_.createDBTable('NOTES',  self._notesdb_sqlcols_)
        else:
            print('INITDB: connecting notesDB')
            self._notesdb_.connectSQLiteDB(self._notesdb_path_)

    def getMetaProteomeData(self, fname):
        assert os.path.exists(fname), "\nERROR ... \n{} does not exist".format(fname)
        F = open(fname, 'r')
        dat = json.load(F)
        F.close()
        mods = {}
        for a in dat:
            if not a.startswith('__'):
                mods[a] = cbmpy.readSBML3FBC(dat[a]['sbml_out'])
                mods[a].setId(cbmpy.CBXML.formatSbmlId(str(a)))
                mods[a].setName(str(a))
                print(mods[a].getId(), mods[a].getName())
        sres = dat['__metaproteome__']['search_results']
        g2react = {}
        for t in sres:
            #print(t, sres[t])
            if sres[t] is not None:
                for g in sres[t]:
                    if g in dat['__idx__']:
                        rids = dat[dat['__idx__'][g]]['gene2reaction'][g]
                        #print(rids)
                        if g in g2react:
                            print('WARNING: {} in G2REACT'.format(g))
                        else:
                            a = [mods[dat['__idx__'][g]].getReaction(r) for r in rids]
                            #print(a)
                            g2react[g] = a
                    else:
                        print('WARNING: {} not in _IDX_: gene defined that does not exist in a GPR association').format(g)
            #print('')

        self._DAT_MODELS = mods
        self._DAT_SEARCH_RES = sres
        self._DAT_G2REACT = g2react
        self._DAT_LINK_DICT_ = dat
        if '__notesdb_key__' in self._DAT_LINK_DICT_['__metaproteome__']:
            print('Using notesdb_key', self._DAT_LINK_DICT_['__metaproteome__']['__notesdb_key__'])
        else:
            self._DAT_LINK_DICT_['__metaproteome__']['__notesdb_key__'] = str(time.time())
            print('Updating results with notesdb_key')
            self.func_saveResultsFile()
        self._DAT_NOTESDB_KEY_ = self._DAT_LINK_DICT_['__metaproteome__']['__notesdb_key__']

    def copyFunc(self):
        return

    def pasteFunc(self):
        return

    def createMenus(self):
        menubar = self.appwindow.menuBar()
        #menubar.setMinimumWidth(300)
        #menubar.setNativeMenuBar(False)

        # menu actions
        # quit
        quitApp = QAction('Quit', self)
        quitApp.setShortcut('Ctrl+Q')
        quitApp.triggered.connect(self.appwindow.close)

        sbmlApp2 = QAction('Export SBML FBCv2', self)
        sbmlApp2.setShortcut('Ctrl+0')
        sbmlApp2.triggered.connect(self.menu_exportSBML2)
        sbmlApp2.setDisabled(True)

        sbmlApp1 = QAction('Export SBML FBCv1', self)
        sbmlApp1.setShortcut('Ctrl+1')
        sbmlApp1.triggered.connect(self.menu_exportSBML1)
        sbmlApp1.setDisabled(True)

        sbmlAppArch = QAction('Export COMBINE Archive', self)
        sbmlAppArch.setShortcut('Ctrl+2')
        sbmlAppArch.triggered.connect(self.menu_exportCOMBINE)
        sbmlAppArch.setDisabled(True)

        sbmlApp0 = QAction('Export COBRA SBML', self)
        sbmlApp0.setShortcut('Ctrl+3')
        sbmlApp0.triggered.connect(self.menu_exportSBML0)
        sbmlApp0.setDisabled(True)

        excelApp = QAction('Export Model as Excel', self)
        excelApp.setShortcut('Ctrl+L')
        excelApp.triggered.connect(self.menu_exportExcel)
        excelApp.setDisabled(True)

        tblxApp = QAction('Export Table as CSV', self)
        tblxApp.setShortcut('Ctrl+T')
        tblxApp.triggered.connect(self.menu_exportTableToCSV)
        tblxApp.setDisabled(True)

        sconApp = QAction('Identify/convert SBML', self)
        sconApp.setShortcut('Ctrl+I')
        sconApp.triggered.connect(self.menu_convertSBML)
        sconApp.setDisabled(False)

        resetApp = QAction('New analysis (reset)', self)
        #resetApp.setShortcut('Ctrl+L')
        #resetApp.isEnabled(False)
        resetApp.triggered.connect(self.menu_resetGUI)
        resetApp.setDisabled(True)

        aboutApp = QAction('About MetaDraft', self)
        aboutApp.triggered.connect(self.menu_helpAbout)

        enableBenchmarkApp = QAction('Use benchmark sequence', self)
        #enableBenchmarkApp.setShortcut('Ctrl+B')
        enableBenchmarkApp.triggered.connect(self.menu_enableBenchmark)

        addSeqPlusModelApp = QAction('Create template model', self)
        #addSeqPlusModelApp.setShortcut('Ctrl+B')
        addSeqPlusModelApp.triggered.connect(self.menu_addSeqPlusModel)

        self.enableOptimizationApp = QAction('Use ID optimization', self, checkable=True)
        self.enableOptimizationApp.setChecked(True)

        self.userMetaDefMenu = QMenu('MetaProteomes', self)
        userMetaDefApp_add = QAction('Add', self)
        self.userMetaDefMenu.addAction(userMetaDefApp_add)

        userMetaDefApp_export = QAction('Export', self)
        self.userMetaDefMenu.addAction(userMetaDefApp_export)

        userMetaDefApp_del = QAction('Delete', self)
        self.userMetaDefMenu.addAction(userMetaDefApp_del)

        self.userMetaDefMenu.addSeparator()
        for m in self._CONFIG_['users'][self.func_getCurrentUser()]['metaproteomes']:
            self.userMetaDefMenu.addAction(QAction(m, self))
        self.userMetaDefMenu.triggered[QAction].connect(self.menu_userMetaDefApp)

        configMenuItem = QAction('Config Homolgy', self)
        configMenuItem.triggered.connect(self.menu_configMenuItem)

        if self.__SYSLOG_ENABLED__:
            viewSyslogApp = QAction('View SysLog', self)
            viewSyslogApp.triggered.connect(self.menu_viewSyslogApp)

        self.savedSessionMenu = QMenu('Saved sessions', self)
        saveSessionApp  = QAction('Save current session', self)
        saveSessionApp.triggered.connect(self.menu_saveSession)

        clearSessionsApp  = QAction('Clear all sessions', self)
        clearSessionsApp.triggered.connect(self.menu_clearSessions)

        menuFile = menubar.addMenu('&File')
        menuFile.addAction(resetApp)
        menuFile.addSeparator()
        menuFile.addAction(sbmlApp2)
        menuFile.addAction(excelApp)
        menuFile.addAction(sbmlAppArch)
        menuFile.addSeparator()
        menuFile.addAction(sbmlApp1)
        menuFile.addAction(sbmlApp0)
        menuFile.addSeparator()
        menuFile.addAction(tblxApp)
        menuFile.addSeparator()
        menuFile.addAction(sconApp)
        menuFile.addSeparator()
        menuFile.addAction(quitApp)

        menuTools = menubar.addMenu('&Build options')
        menuTools.addAction(enableBenchmarkApp)
        menuTools.addAction(addSeqPlusModelApp)
        #menuTools.addSeparator()
        menuTools.addAction(self.enableOptimizationApp)
        menuTools.addSeparator()
        menuTools.addMenu(self.userMetaDefMenu)
        menuTools.addSeparator()
        menuTools.addAction(configMenuItem)


        self.menuToolsM = menubar.addMenu('&Model options')

        model_actions = ['Summary report', 'Gene report', 'Reaction report', 'Metabolite report', 'Export unmatched genes (FASTA)', 'Export unselected genes (FASTA)', 'Export model notes (CSV)']
        for m in model_actions:
            self.menuToolsM.addAction(QAction(m, self))
        self.menuToolsM.triggered[QAction].connect(self.menu_modelTools)
        self.menuToolsM.setDisabled(True)

        self.menuSessions = menubar.addMenu('&Sessions')
        self.menuSessions.setDisabled(True)
        self.menuSessions.addAction(saveSessionApp)
        self.menuSessions.addMenu(self.savedSessionMenu)
        self.menuSessions.addAction(clearSessionsApp)
        self.menuTools = menuTools

        self.sbmlApp0 = sbmlApp0
        self.sbmlApp1 = sbmlApp1
        self.sbmlApp2 = sbmlApp2
        self.sbmlAppArch = sbmlAppArch
        self.excelApp = excelApp
        self.resetApp = resetApp
        self.tblxApp = tblxApp

        menuHelp = menubar.addMenu('&Help')
        if self.__SYSLOG_ENABLED__:
            menuHelp.addAction(viewSyslogApp)
        menuHelp.addAction(aboutApp)

    def menu_savedSessionMenu(self):
        cstate = self._DAT_LINK_DICT_['__metaproteome__']['selection_state']
        cstate_keys = list(cstate.keys())
        cstate_keys.sort(reverse=True)
        if '__tmp__' in cstate_keys:
            cstate_keys.remove('__tmp__')
        self.savedSessionMenu.clear()

        self.sessionSignalMapper = QtCore.QSignalMapper(self)
        self._session_actions_ = []
        for k in range(len(cstate_keys)):
            self._session_actions_.append(QAction(cstate_keys[k], self))
            self.sessionSignalMapper.setMapping(self._session_actions_[k], cstate_keys[k])
            self._session_actions_[k].triggered.connect(self.sessionSignalMapper.map)
            self.savedSessionMenu.addAction(self._session_actions_[k])
        self.sessionSignalMapper.mapped[str].connect(self.menu_loadSession)
        if len(cstate_keys) > 0:
            self.status_bar.showMessage('Session saved as: {}.'.format(cstate_keys[0]))

    @pyqtSlot(str)
    def menu_loadSession(self, state):
        #print(state)
        self.func_loadSelectionState(state)
        self.func_setSelectionState()

    @pyqtSlot()
    def menu_clearSessions(self):
        self.savedSessionMenu.clear()
        cstate = self._DAT_LINK_DICT_['__metaproteome__']['selection_state']
        for k in list(cstate.keys()):
            if k != '__tmp__':
                cstate.pop(k)
        self.func_saveSelectionState(True, False, False)
        self.status_bar.showMessage('Saved sessions cleared.')

    @pyqtSlot()
    def menu_saveSession(self):
        self.func_saveSelectionState(True, True, True)
        self.menu_savedSessionMenu()

    def menu_buildAll(self):
        self._updateGeneMap_()
        self._update_Reactions_()
        self._update_Metabolites_()

    @pyqtSlot()
    def menu_helpAbout(self):

        title = "About MetaDraft"
        msg = "This is the MetaDraft version {}-({}), https://github.com/SystemsBioinformatics/metadraft. MetaDraft uses PySCeS-CBMPy (http://cbmpy.sourceforge.net) technology, both part of the MetaToolkit project.\n\n".format(__version__, cbmpy.__version__)
        msg += "(c) Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, 2015-2017. All rights reserved\n"
        msg += "\nFor help and support please contact Brett Olivier:\nemail: b.g.olivier@vu.nl\n"
        if HAVE_QT4:
            qtv = 'Qt4'
        else:
            qtv = 'Qt5'
        msg += "\nYou are using Py{}.\n".format(qtv)
        self.widgetMsgBox(QMessageBox.Information, title, msg)

    @pyqtSlot(QAction)
    def menu_modelTools(self, q):
        rpt = str(q.text())
        print(rpt)
        if rpt == 'Summary report':
            self.func_generateSummaryReport()
            self.widget_displayReport(self.func_formatSummaryReport())
        elif rpt == 'Gene report':
            self.func_generateGeneReport()
            self.widget_displayReport(self.func_formatGeneReport())
        elif rpt == 'Reaction report':
            self.func_generateReactionReport()
            self.widget_displayReport(self.func_formatReactionReport())
        elif rpt == 'Metabolite report':
            self.func_generateMetaboliteReport()
            self.widget_displayReport(self.func_formatMetaboliteReport())
        elif rpt == 'Full report':
            self.func_generateGeneReport()
            self.func_generateReactionReport()
            self.func_generateMetaboliteReport()
        elif rpt == 'Export unmatched genes (FASTA)':
            self.func_exportUnmatchedGenes()
        elif rpt == 'Export unselected genes (FASTA)':
            self.func_exportUnselectedGenes()
        elif rpt == 'Export model notes (CSV)':
            self.func_exportNotesDB()
        else:
            print('NOTHING TO REPORT SIRE!')

    def func_exportNotesDB(self):
        filename = str(self.saveFile('Export notes database', self._history_save_dir_, '*.csv'))
        try:
            #self._notesdb_.dumpTableToTxt('NOTES', filename)
            self._notesdb_.dumpTableToCSV('NOTES', filename)
        except IOError:
            print('exportNotesDB, no filename')

    def func_exportUnmatchedGenes(self):
        LD = self._DAT_LINK_DICT_['__metaproteome__']
        unmatched = []
        for g in LD['search_results']:
            if LD['search_results'][g] is None:
                unmatched.append(g[len(self.gene_prefix):])
        #print(unmatched)
        fname = self.openFile('Load Input File', self._history_open_dir_, 'Supported (*.gbk *.gbff *.gb *.fasta *.faa *.fa);;FASTA (*.fasta *.faa *.fa);;GenBank (*.gbk *.gbff *.gb)')
        seq = []
        seqout = {}
        if fname.endswith('.gbk') or fname.endswith('.gb') or fname.endswith('.gbff'):
            GBFile = open(fname, 'r')
            GBcds = biotools.Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
            for cds in GBcds:
                if cds.seq != None:
                    cds.id = cds.name
                    #cds.description = ''
                    #if gene_prefix is not None:
                        #cds.id = gene_prefix + cds.id
                    seqout[cds.name] = cds
            GBFile.close()
            #seq = biotools.SeqIO.read(fname, 'genbank')

            #for fe_ in seq.features:
                #if fe_.type == 'CDS':
                    #if fe_.qualifiers['locus_tag'][0] in unmatched:
                        #prt = fe_.qualifiers['translation'][0]
                        #gid = fe_.qualifiers['locus_tag'][0]
                        #prt = biotools.Bio.SeqRecord.SeqRecord(biotools.Bio.Seq.Seq(prt,\
                                                                                    #biotools.Bio.Alphabet.ProteinAlphabet()),\
                                                               #id=gid, name=gid, description='')
                        #seqout[gid] = prt
        else:
            try:
                seq = biotools.SeqIO.parse(fname, 'fasta')
                for s in seq:
                    if s.id in unmatched:
                        seqout[s.id] = s
            except:
                print('FASTA file read error')
                self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fname)))
        #print(seqout)
        #print(len(seqout))
        biotools.writeFASTA(fname+'.unmatched.fasta', seqout, paranoid_style=False)

    def func_exportUnselectedGenes(self):
        LD = self._DAT_LINK_DICT_['__metaproteome__']

        self._updateGeneMap_()
        unmatched = []
        for g in LD['search_results']:
            if LD['search_results'][g] is None:
                unmatched.append(str(g[len(self.gene_prefix):]))

        unselected = []
        for g in self._gene_selected_map_:
            if not self._gene_selected_map_[g]:
                g1 = g.split(self.id_sep)[0][len(self.gene_prefix):]
                if g1 not in unmatched:
                    unselected.append(g1)

        fname = self.openFile('Load Input File', self._history_open_dir_, 'Supported (*.gbk *.gbff *.gb *.fasta *.faa *.fa);;FASTA (*.fasta *.faa *.fa);;GenBank (*.gbk *.gbff *.gb)')
        seq = []
        seqout = {}
        if fname.endswith('.gbk') or fname.endswith('.gb') or fname.endswith('.gbff'):
            GBFile = file(fname, 'r')
            GBcds = biotools.Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
            for cds in GBcds:
                if cds.seq != None:
                    cds.id = cds.name
                    seqout[cds.name] = cds
            GBFile.close()
        else:
            try:
                seq = biotools.SeqIO.parse(fname, 'fasta')
                for s in seq:
                    if s.id in unselected:
                        seqout[s.id] = s
            except:
                print('FASTA file read error')
                self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fname)))
        biotools.writeFASTA(fname+'.unselected.fasta', seqout, paranoid_style=False)

    def func_generateSummaryReport(self):
        self.func_generateGeneReport()
        self.func_generateReactionReport()
        self.func_generateMetaboliteReport()

    def func_generateGeneReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']
        LD['reports']['genes']['unmatched'] = []
        LD['reports']['genes']['selected'] = []
        LD['reports']['genes']['unselected'] = []
        for g in LD['search_results']:
            if LD['search_results'][g] is None:
                LD['reports']['genes']['unmatched'].append(g[len(self.gene_prefix):])
        for gidx in range(self.table_gene.rowCount()):
            igene = str(self.table_gene.item(gidx, 0).text()).strip()[len(self.gene_prefix):]
            mgene = str(self.table_gene.item(gidx, 1).text()).strip()
            score = str(self.table_gene.item(gidx, 2).text()).strip()
            org = str(self.table_gene.item(gidx, 4).text()).strip()
            gdat = (igene, mgene, score, org)
            if igene not in LD['reports']['genes']['unmatched']:
                if self.table_gene.item(gidx, 3).checkState():
                    LD['reports']['genes']['selected'].append(gdat)
                else:
                    LD['reports']['genes']['unselected'].append(gdat)
        QApplication.restoreOverrideCursor()


    def func_generateReactionReport(self): #r, n, o, s
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']
        LD['reports']['reactions']['selected'] = []
        LD['reports']['reactions']['unselected'] = []
        for ridx in range(self.table_reaction.rowCount()):
            react = str(self.table_reaction.item(ridx, 0).text()).strip()
            try:
                name = str(self.table_reaction.item(ridx, 1).text()).strip()
            except UnicodeEncodeError:
                name = 'invalid name'
            org = str(self.table_reaction.item(ridx, 3).text()).strip()
            source = str(self.table_reaction.item(ridx, 5).text()).strip()
            rdat = (react, name, org, source)
            if self.table_reaction.item(ridx, 2).checkState():
                LD['reports']['reactions']['selected'].append(rdat)
            else:
                LD['reports']['reactions']['unselected'].append(rdat)
        QApplication.restoreOverrideCursor()


    def func_generateMetaboliteReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']
        LD['reports']['metabolites']['selected'] = []
        LD['reports']['metabolites']['unselected'] = []
        for midx in range(self.table_metab.rowCount()):
            metab = str(self.table_metab.item(midx, 0).text()).strip()
            try:
                name = str(self.table_metab.item(midx, 1).text()).strip()
            except UnicodeEncodeError:
                name = 'invalid name'
            org = str(self.table_metab.item(midx, 2).text()).strip()
            mdat = (metab, name, org)
            LD['reports']['metabolites']['selected'].append(mdat)
            #if self.table_metab.item(midx, 3).checkState():
                #LD['reports']['metabolites']['selected'].append(mdat)
            #else:
                #LD['reports']['metabolites']['unselected'].append(mdat)
        QApplication.restoreOverrideCursor()

    def func_formatSummaryReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']['reports']
        gene_stats = {'org' : {},
                      'score' : {}
                     }
        # summarize gene info
        total_sel = 0
        for g, m, s, o in LD['genes']['selected']:
            total_sel += 1
            if o not in gene_stats['org']:
                gene_stats['org'][o] = 1
                gene_stats['score'][o] = float(s)
            else:
                gene_stats['org'][o] += 1
                gene_stats['score'][o] += float(s)

        #print(gene_stats)

        gene_tbl = []
        for o in gene_stats['org']:
            gene_tbl.append((gene_stats['org'][o], o))
        gene_tbl.sort()
        gene_tbl.reverse()

        # summarize reaction info
        total_sel_react = 0
        react_stats = {'org' : {}}
        for r, n, o, s in LD['reactions']['selected']:
            total_sel_react += 1
            if o not in react_stats['org']:
                react_stats['org'][o] = 1
            else:
                react_stats['org'][o] += 1

        #print(react_stats)

        react_tbl = []
        for o in react_stats['org']:
            react_tbl.append((react_stats['org'][o], o))
        react_tbl.sort()
        react_tbl.reverse()

        # summarize metabolite info
        total_sel_metab = 0
        metab_stats = {'org' : {}}
        for m, n, o in LD['metabolites']['selected']:
            total_sel_metab += 1
            if o not in metab_stats['org']:
                metab_stats['org'][o] = 1
            else:
                metab_stats['org'][o] += 1

        #print(metab_stats)

        metab_tbl = []
        for o in metab_stats['org']:
            metab_tbl.append((metab_stats['org'][o], o))
        metab_tbl.sort()
        metab_tbl.reverse()

        html = """
        <html>
        <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>MetaDraft Summary Report</title>
        </head>
        <body>
            <h2 id="page-top">MetaDraft Summary Report</h2>
            <p><strong>Input file:</strong> {}<br/><strong>Metaproteome:</strong> {}<br/><strong>Result file:</strong> {}<br/><strong>Report date:</strong> {}<br/><strong>MetaDraft version:</strong> {}</p>
            <h3>Gene statistics</h3>
            <table width="80%" border="0" cellspacing="1" cellpadding="1">
                <tbody>
        """.format(os.path.split(self._DAT_LINK_DICT_['__metaproteome__']['input_fasta'])[-1].replace('.in.','.'),\
                                          self._DAT_LINK_DICT_['__metaproteome__']['__fullname__'],\
                                          self.result_file, time.strftime('%Y-%m-%d %H:%M'), __version__)

        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Model', 'Percentage', 'Genes', 'Avg. score')
        for o in gene_tbl:
            html += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o[1], round(float(o[0])/total_sel*100.0, 2), o[0],\
                                                                        round(float(gene_stats['score'][o[1]])/gene_stats['org'][o[1]], 3))
        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Total', '', total_sel , '')
        html += """
            </tbody>
        </table>"""

        html += """    <h3>Reaction statistics</h3>
            <table width="80%" border="0" cellspacing="1" cellpadding="1">
                <tbody>
        """

        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Model', 'Percentage', 'Reactions', '')
        for o in react_tbl:
            html += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o[1], round(float(o[0])/float(total_sel_react)*100.0, 2), o[0], '')
        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Total', '', total_sel_react , '')

        html += """
            </tbody>
        </table>"""

        html += """    <h3>Metabolite statistics</h3>
            <table width="80%" border="0" cellspacing="1" cellpadding="1">
                <tbody>
        """

        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Model', 'Percentage', 'Metabolites', '')
        for o in metab_tbl:
            html += '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o[1], round(float(o[0])/float(total_sel_metab)*100.0, 2), o[0], '')
        html += '<tr><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td><td><strong>{}</strong></td></tr>'.format('Total', '', total_sel_metab , '')

        html += """
            </tbody>
        </table>"""

        html += """
        {}
        </body>
        </html>
        """
        cp = self.metadraft_rpt_footer
        QApplication.restoreOverrideCursor()
        return html.format(cp)

    def func_formatGeneReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']['reports']
        nmg = ', '.join(LD['genes']['unmatched'])
        selg = '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Source', 'Target', 'Score', 'Model')
        def buildAnnotation(m, g, o):
            annot = self.buildHtmlStringsGene(m, g)
            annot = annot.replace('<table cellpadding=', '<br/><table id="{}_{}_{}" width="80%" cellpadding='.format(o, g, m))
            annot = annot.replace('</caption>',' (<a href="#page-top">top</a>)</caption>')
            annot = annot.replace('</table>','</table>\n')
            return annot
        mgene = []
        annotG = ''
        for g, m, s, o in LD['genes']['selected']:
            if g not in mgene:
                mgene.append(g)
            selg += '<tr><td><a href="#{}_{}_{}">{}</a></td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o, g, m, g, m, s, o)
            annotG += buildAnnotation(m, g, o)
        nselg = '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Source', 'Target', 'Score', 'Model')
        tgene = []
        for g, m, s, o in LD['genes']['unselected']:
            if g not in tgene and g not in mgene:
                tgene.append(g)
            nselg += '<tr><td><a href="#{}_{}_{}">{}</a></td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o, g, m, g, m, s, o)
            annotG += buildAnnotation(m, g, o)
        #print(mgene, len(mgene))
        cp = self.metadraft_rpt_footer
        QApplication.restoreOverrideCursor()
        return metadraftreports.gene_report_template(__version__, '{}:{}'.format(len(mgene), len(tgene+mgene)), os.path.split(self._DAT_LINK_DICT_['__metaproteome__']['input_fasta'])[-1].replace('.in.','.'),\
                                          self._DAT_LINK_DICT_['__metaproteome__']['__fullname__'], nmg, selg, nselg, annotG, cp)


    def func_formatReactionReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']['reports']
        annotR = ''
        def buildAnnotation(r, o):
            R = self._DAT_MODELS[o].getReaction(r)
            annot = self.buildHtmlStringsReaction(R)
            annot = annot.replace('<table>','<hr/><table id="{}_{}"><caption>{} (<a href="#page-top">top</a>)</caption>'.format(o, r, r))
            annot = annot.replace('<html><head></head><body>','')
            annot = annot.replace('</body></html>','\n')
            miriam = R.getMIRIAMannotations()
            r_html = ''
            if miriam != None:
                r_html += "<tr><td colspan=\"2\" align=\"left\">RDF references (opens in browser)</td></tr>"
                for m in miriam:
                    if len(miriam[m]) > 0:
                        for u in range(len(miriam[m])):
                            r_html += "<tr><td>{}</td><td><a href=\"{}\">{}</a></td></tr>".format(m, miriam[m][u], miriam[m][u])
            annot = annot.replace('</table>', '\n{}\n</table>'.format(r_html))
            return annot

        selg = '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Reaction', 'Name', 'Model', 'Source')
        nreact = 0
        for r, n, o, s in LD['reactions']['selected']:
            nreact += 1
            selg += '<tr><td><a href="#{}_{}">{}</a></td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o, r, r, n, o, s)
            annotR += buildAnnotation(r, o)
        nselg = '<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Reaction', 'Name', 'Model', 'Source')
        for r, n, o, s in LD['reactions']['unselected']:
            nselg += '<tr><td><a href="#{}_{}">{}</a></td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(o, r, r, n, o, s)
            annotR += buildAnnotation(r, o)
        cp = self.metadraft_rpt_footer
        QApplication.restoreOverrideCursor()
        return metadraftreports.reaction_report_template(__version__, nreact, os.path.split(self._DAT_LINK_DICT_['__metaproteome__']['input_fasta'])[-1].replace('.in.','.'),\
                                          self._DAT_LINK_DICT_['__metaproteome__']['__fullname__'], selg, nselg, annotR, cp)


    def func_formatMetaboliteReport(self):
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        LD = self._DAT_LINK_DICT_['__metaproteome__']['reports']
        selg = '<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Metabolite', 'Name', 'Model')
        def buildAnnotation(m, o):
            annot = self.buildHtmlStringsMetab(m)
            annot = annot.replace('<table cellpadding="5"','<hr/><table cellpadding="5" width="80%" id="{}_{}"'.format(o, m))
            annot = annot.replace('<caption></caption>','<caption>{} (<a href="#page-top">top</a>)</caption>'.format(m))
            annot = annot.replace('<html><body>','')
            annot = annot.replace('</body></html>','\n')
            return annot
        annotM = ''
        nmetab = 0
        for m, n, o in LD['metabolites']['selected']:
            nmetab += 1
            selg += '<tr><td><a href="#{}_{}">{}</a></td><td>{}</td><td>{}</td></tr>'.format(o, m, m, n, o)
            annotM += buildAnnotation(m, o)
        nselg = '<br/>'
        #nselg = '<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format('Metabolite', 'Name', 'Model')
        #for m, n, o in LD['metabolites']['unselected']:
            #nselg += '<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format(m, n, o)
        cp = self.metadraft_rpt_footer
        QApplication.restoreOverrideCursor()
        return metadraftreports.metabolite_report_template(__version__, nmetab, os.path.split(self._DAT_LINK_DICT_['__metaproteome__']['input_fasta'])[-1].replace('.in.','.'),\
                                          self._DAT_LINK_DICT_['__metaproteome__']['__fullname__'], selg, nselg, annotM, cp)

    def widget_displayReport(self, html):
        try:
            html = html.encode('utf-8', 'ignore').strip()
        except:
            try:
                html = html.decode('utf-8', 'ignore')
            except Exception as ex:
                print('Could not encode report')
                print(ex)

        self.widget_reportViewApp = QWidget()
        self.widget_reportViewApp.setWindowTitle('MetaDraft Report')
        ppos = self.func_getNewPopupWindowCoords()
        self.widget_reportViewApp.setGeometry(QtCore.QRect(ppos[0], ppos[1], 600, 300))
        self.widget_reportViewApp.setWindowModality(Qt.ApplicationModal)

        textWindow = QTextBrowser(self.widget_reportViewApp)
        textWindow.setOpenLinks(True)
        textWindow.setReadOnly(True)
        try:
            textWindow.setText(html)
        except TypeError:
            textWindow.setText(str(html))

        @pyqtSlot()
        def exportHTML():
            #html = textWindow.toHtml()
            #html = html.encode('utf-8', 'ignore').strip()
            filename = str(self.saveFile('Save report', self._history_save_dir_, '*.html'))
            try:
                F = open(filename, 'w')
                F.write(html)
                F.close()
            except IOError:
                print('exportHTML, no filename')
            self.widget_reportViewApp.close()

        @pyqtSlot()
        def openBrowser():
            path = os.path.abspath(os.path.join(self._tmpDir_, '_view_rpt_.html'))
            url = 'file://' + path
            with open(path, 'w') as f:
                try:
                    f.write(html)
                except TypeError:
                    f.write(str(html))
            webbrowser.open_new_tab(url)
            self.widget_reportViewApp.close()

        savebut = QPushButton(self.widget_reportViewApp)
        savebut.setText('Export (HTML)')
        savebut.clicked.connect(exportHTML)

        browsebut = QPushButton(self.widget_reportViewApp)
        browsebut.setText('View in browser')
        browsebut.clicked.connect(openBrowser)

        #def addNotes():
            #<!--<h3>Notes</h3>
            #<p></p>-->
            #<h3>Genes</h3>
            #pass
        #addnotes = QPushButton(self.widget_reportViewApp)
        #addnotes.setText('Add notes')
        #addnotes.clicked.connect(addNotes)


        layout = QGridLayout(self.widget_reportViewApp)
        layout.setSpacing(10)
        layout.addWidget(textWindow, 0, 0, 4, 2)
        layout.addWidget(browsebut, 4, 0, 1, 1)
        layout.addWidget(savebut, 4, 1, 1, 1)

        self.widget_reportViewApp.show()

    @pyqtSlot()
    def menu_enableBenchmark(self):
        if not self.bp_btn_out.isEnabled():
            self.bp_text_out.setEnabled(True)
            self.bp_label_out.setEnabled(True)
            self.bp_btn_out.setEnabled(True)
            self.bp_targetOpen2()
        else:
            self.bp_text_out.setText('')
            self.bp_text_out.setDisabled(True)
            self.bp_label_out.setDisabled(True)
            self.bp_btn_out.setDisabled(True)

    @pyqtSlot()
    def menu_configMenuItem(self):
        print(bionoid.CONFIGKEYS)
        self.widget_config = ConfigPanelWidgetINP(bionoid, 'CONFIGKEYS')
        print(bionoid.CONFIGKEYS)




    @pyqtSlot(QAction)
    def menu_userMetaDefApp(self, q):
        mset = str(q.text())
        if mset == 'Add':
            self.menu_userMetaDefApp_add()
        elif mset == 'Delete':
            self.menu_userMetaDefApp_del()
        elif mset == 'Export':
            self.menu_userMetaDefApp_export()
        else:
            metap = self._CONFIG_['users'][self.func_getCurrentUser()]['metaproteomes'][mset]
            selected = metap.split('&nbsp;')

            all_items = [str(self.bp_lview.item(i).text()) for i in range(self.bp_lview.count())]
            selected = [i for i in selected if i in all_items]
            print(selected)

            pre_selected = len(selected)
            self.bp_lview.selectAll()
            list_items = self.bp_lview.selectedItems()
            for i in list_items:
                if str(i.text()) not in selected:
                    selected.append(str(i.text()))
            self.bp_lview.clear()
            self.bp_lview.addItems(selected)
            for i in range(pre_selected):
                self.bp_lview.item(i).setSelected(True)
            self.bp_lview.update()
            self.bp_btn_blast.setEnabled(True)
            self.bp_btn_blast.update()

    def menu_userMetaDefApp_add(self):
        s_items = self.bp_lview.selectedItems()
        if len(s_items) > 0:
            s_items = [str(a.text()) for a in s_items]
            s_items = '&nbsp;'.join(s_items)
            mset, ok = QInputDialog.getText(self, 'Input Dialog',
                                                  'Enter MetaProteome name:')
            if ok:
                self._CONFIG_['users'][self.func_getCurrentUser()]['metaproteomes'][str(mset)] = s_items
                self.userMetaDefMenu.addAction(QAction(str(mset), self))
                self._writeConfig()

    def func_getNewPopupWindowCoords(self):
        ppos = self.mapToGlobal(self.pos())
        size = self.size()
        w, h = int(size.width()/2 - 200), int(size.height()/2 - 100)
        x, y, = int(ppos.x()) + w, int(ppos.y()) + h
        return x, y

    def menu_userMetaDefApp_del(self):
        print('userMetaDefApp_del')
        self.widget_userMetaDefApp_del = QWidget()
        self.widget_userMetaDefApp_del.setWindowTitle('Delete user defined list(s)')
        ppos = self.func_getNewPopupWindowCoords()
        self.widget_userMetaDefApp_del.setGeometry(QtCore.QRect(ppos[0], ppos[1], 400, 200))

        list_widge = QListWidget(parent=self.widget_userMetaDefApp_del)
        #list_widge.setDragDropMode(QAbstractItemView.InternalMove)
        list_widge.setSelectionMode(QAbstractItemView.ExtendedSelection)
        litems = {}
        for action in self.userMetaDefMenu.actions():
            txt = str(action.text())
            if txt not in ['Add', 'Export', 'Delete', '']:
                litems[txt] = action
        list_widge.addItems(litems.keys())

        @pyqtSlot()
        def exit():
            self.widget_userMetaDefApp_del.close()

        @pyqtSlot()
        def delete():
            delist = list_widge.selectedItems()
            if len(delist) > 0:
                for i in delist:
                    i = str(i.text())
                    self.userMetaDefMenu.removeAction(litems[i])
                    self._CONFIG_['users'][self.func_getCurrentUser()]['metaproteomes'].pop(i)
                self._writeConfig()
            self.widget_userMetaDefApp_del.close()

        layout = QGridLayout(self.widget_userMetaDefApp_del)
        layout.setSpacing(10)

        delbut = QPushButton(self.widget_userMetaDefApp_del)
        delbut.setText('Delete and close')
        delbut.clicked.connect(delete)

        exbut = QPushButton(self.widget_userMetaDefApp_del)
        exbut.setText('Close')
        exbut.clicked.connect(exit)

        layout.addWidget(list_widge, 0, 0, 3, 2)
        layout.addWidget(delbut, 3, 0)
        layout.addWidget(exbut, 3, 1)
        self.widget_userMetaDefApp_del.setWindowModality(Qt.ApplicationModal)
        self.widget_userMetaDefApp_del.show()


    @pyqtSlot(QtCore.QPoint)
    def bp_lviewRightClicked(self, QPos):

        def deleteModel():
            print(item)
            path = os.path.join(self.seqplus_files, '{}.seqplus.xml'.format(item[0]))
            print(path)
            if os.path.exists(path):
                reply = QMessageBox.question(self, 'Message',\
                                                   "Are you sure you want to delete:\n{}?".format(item[0]),\
                                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if reply == QMessageBox.Yes:
                    os.remove(path)
                    os.remove(path.replace('.xml','.json'))
                    print('Deleted files: \"{}\", \"{}\"'.format(path, path.replace('.xml','.json')))
                    self.bp_lview.clear()
                    datF = []
                    for f in os.listdir(self.seqplus_files):
                        if f.endswith('.seqplus.xml'):
                            datF.append(f.replace('.seqplus.xml',''))
                    datF.sort()
                    self.bp_lview.addItems(datF)

        def exportNonGPR():
            print('export nonGPR')
            print(item)
            QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
            exm = cbmpy.CBModel.Model('metadraft_nongpr_reactions')
            notes = '<p>Metadraft nonGPR reactions for models:</p><p><br/>'
            idGlob = {}
            for m in item:
                reactions = []
                notes += '{}.seqplus.xml<br/>'.format(m)
                modelF = os.path.join(self.seqplus_files, m + '.seqplus.xml')
                cmod = cbmpy.readSBML3FBC(modelF)
                reactions = self.func_findNonGprReactionsForModel(cmod)
                for rid in reactions:
                    if rid not in idGlob:
                        idGlob[rid] = 0
                        cbmpy.CBMultiModel.copyReaction(cmod, exm, rid)
                    else:
                        idGlob[rid] = idGlob[rid] + 1
                        altrid='{}_{}'.format(rid, idGlob[rid])
                        idGlob[altrid] = idGlob[rid]
                        #print(modelF, rid, altrid)
                        cbmpy.CBMultiModel.copyReaction(cmod, exm, rid, altrid)
            notes += '</p>'
            print(len(reactions))
            exm.setNotes(notes)
            for s in exm.species:
                s.unsetBoundary()
                s.setCompartmentId('cell')
            QApplication.restoreOverrideCursor()
            filename = str(self.saveFile('Export NonGPR SBML', self.sbml_save_dir, '*.xml'))
            cbmpy.writeSBML3FBCV2(exm, filename)

        item = self.bp_lview.selectedItems()
        self.widget_lview_rclickmenu = QMenu()
        if len(item) > 1:
            item = [str(i.text()) for i in item]
            #self.widgetMsgBox(QMessageBox.Information, 'Template model menu', 'Please select only one item to enable right click menu.')
            menu_item2 = self.widget_lview_rclickmenu.addAction("Export non-gpr reac.")
            menu_item2.triggered.connect(exportNonGPR)
        else:
            item = [str(item[0].text())]
            menu_item2 = self.widget_lview_rclickmenu.addAction("Export non-gpr reac.")
            menu_item2.triggered.connect(exportNonGPR)
            self.widget_lview_rclickmenu.addSeparator()
            menu_item = self.widget_lview_rclickmenu.addAction("Delete model")
            menu_item.triggered.connect(deleteModel)

        parentPosition = self.bp_lview.viewport().mapToGlobal(QPos)
        self.widget_lview_rclickmenu.move(parentPosition)
        self.widget_lview_rclickmenu.show()

    def menu_userMetaDefApp_export(self):
        print('userMetaDefApp_export')
        items = self.bp_lview.selectedItems()
        items2 = []
        for i in items:
            sp_str = str(i.text()).split(')-(')
            items2.append('{}'.format(sp_str[0][1:]))

        self.optimized_metaproteome = self.enableOptimizationApp.isChecked()
        print('Optimization:', self.optimized_metaproteome)
        self.buildMetaProteomeFromSeqplus(items2, self.seqplus_files,
                                          self.optimized_metaproteome)

        out_file = self.metaproteome_file
        print(out_file)
        #self.metaproteome_file = None
        filename = str(self.saveFile('Export Metaproteome to:', self._history_save_dir_, 'All files (*)'))
        filename = os.path.abspath(filename)
        shutil.move(out_file, filename)
        os.remove(out_file.replace('_metaproteome.fasta','_metalink.json'))
        print(filename)

        #s_items = self.bp_lview.selectedItems()
        #if len(s_items) > 0:
            #s_items = [str(a.text()) for a in s_items]
            #s_items = '&nbsp;'.join(s_items)
            #mset, ok = QInputDialog.getText(self, 'Input Dialog',
                                                  #'Enter MetaProteome name:')
            #if ok:
                #self._CONFIG_['users'][self.func_getCurrentUser()]['metaproteomes'][str(mset)] = s_items
                #self.userMetaDefMenu.addAction(QAction(str(mset), self))
                #self._writeConfig()

    @pyqtSlot()
    def menu_viewSyslogApp(self):
        self.widgetLogView = QTextEdit()
        ppos = self.func_getNewPopupWindowCoords()
        self.widgetLogView.setGeometry(QtCore.QRect(ppos[0], ppos[1], 600, 300))
        self.widgetLogView.setWindowModality(Qt.ApplicationModal)
        self.widgetLogView.setReadOnly(True)
        self.widgetLogView.setWindowTitle('SysLog Viewer ({})'.format(self.log_syslog))
        if self.log_syslog is not None:
            F = open(self.log_syslog, 'r')
            self.widgetLogView.setPlainText(F.read())
            F.close()
            self.widgetLogView.show()

    @pyqtSlot()
    def menu_addSeqPlusModel(self):
        self.widgetAddSeQPlus = QWidget()
        self.widgetAddSeQPlus.setWindowTitle('Create SeQPlus Template Model')
        ppos = self.func_getNewPopupWindowCoords()
        self.widgetAddSeQPlus.setGeometry(QtCore.QRect(ppos[0], ppos[1], 700, 100))

        layout = QGridLayout(self.widgetAddSeQPlus)
        layout.setSpacing(10)

        label_sbml = QLabel(self.widgetAddSeQPlus)
        label_sbml.setText('SBML file:')
        load_sbml = QPushButton(self.widgetAddSeQPlus)
        load_sbml.setText('Open')
        load_sbml.clicked.connect(self.bp_createSeqPlusLoadFileSBML)

        load_gb = QPushButton(self.widgetAddSeQPlus)
        load_gb.setText('Open')
        load_gb.clicked.connect(self.bp_createSeqPlusLoadFilesGB)

        self.widget_seqplus_sbml = QLineEdit(self.widgetAddSeQPlus)
        self.widget_seqplus_sbml.setReadOnly(True)
        self.widget_seqplus_sbml.setText(str(None))

        label_gbank = QLabel(self.widgetAddSeQPlus)
        label_gbank.setText('GenBank file(s):')
        exit = QPushButton(self.widgetAddSeQPlus)
        exit.setText('Exit')
        exit.clicked.connect(self.widgetAddSeQPlus.close)
        self.widget_seqplus_gbank = QLineEdit(self.widgetAddSeQPlus)
        self.widget_seqplus_gbank.setReadOnly(True)
        self.widget_seqplus_gbank.setText(str(None))
        label_key = QLabel(self.widgetAddSeQPlus)
        label_key.setText('Short key with no spaces or dashes, e.g. eco1')
        self.widget_seqplus_key = QLineEdit(self.widgetAddSeQPlus)
        self.widget_seqplus_key.setText(str(None))

        process = QPushButton(self.widgetAddSeQPlus)
        process.setText('Create Model')
        process.clicked.connect(self.bp_createSeqPlusProcess)

        layout.addWidget(label_key, 0, 0)
        layout.addWidget(self.widget_seqplus_key, 0, 1)
        layout.addWidget(label_sbml, 1, 0)
        layout.addWidget(self.widget_seqplus_sbml, 1, 1)
        layout.addWidget(load_sbml, 1, 2)
        layout.addWidget(label_gbank, 2, 0)
        layout.addWidget(self.widget_seqplus_gbank, 2, 1)
        layout.addWidget(load_gb, 2, 2)
        layout.addWidget(exit, 3, 2)
        layout.addWidget(process, 3, 1)
        self.widgetAddSeQPlus.setLayout(layout)

        self.widgetAddSeQPlus.setWindowModality(Qt.ApplicationModal)
        self.widgetAddSeQPlus.show()

    @pyqtSlot()
    def bp_createSeqPlusLoadFileSBML(self):
        fname = str(self.openFile('Load file', self._history_open_dir_, filterex='Supported (*.xml *.sbml);;SBML (*.xml *.sbml)'))
        if fname.endswith('.xml'):
            self.widget_seqplus_sbml.setText(str(fname))
        else:
            self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fname)))
            self.widgetAddSeQPlus.activateWindow()
            return
        self.widget_seqplus_sbml.update()
        self.widgetAddSeQPlus.activateWindow()

    @pyqtSlot()
    def bp_createSeqPlusLoadFilesGB(self):
        fnames = self.openFiles('Load file', self._history_open_dir_, filterex='Supported (*.gbk *.gbff *.gb);;GenBank (*.gbk *.gbff *.gb)')
        outlist = []
        for f in fnames:
            if f.endswith('.gbk') or f.endswith('.gb') or f.endswith('.gbff'):
                outlist.append(f)
        if len(outlist) == 0:
            self.widget_seqplus_gbank.setText('')
            self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fnames)))
            self.widgetAddSeQPlus.activateWindow()
            return
        if len(outlist) == 1:
            self.widget_seqplus_gbank.setText(outlist[0])
        elif len(outlist) > 1:
            self.widget_seqplus_gbank.setText(','.join(outlist))
        self.widget_seqplus_gbank.update()
        self.widgetAddSeQPlus.activateWindow()


    @pyqtSlot()
    def bp_createSeqPlusProcess(self):
        gbank = str(self.widget_seqplus_gbank.text())
        sbml = str(self.widget_seqplus_sbml.text())
        key = str(self.widget_seqplus_key.text())
        ret = 'None'
        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        if gbank == 'None' or sbml == 'None' or key == 'None' or '-' in key:
            QApplication.restoreOverrideCursor()
            title = 'Create SeQPlus Template Model: error'
            msg = 'Model creation error: missing information or  \"-\" in key'
            self.widgetMsgBox(QMessageBox.Critical, title, msg)
            return
        else:
            key = '{}'.format(key)
            sbml_dir, sbml = os.path.split(sbml)[0], os.path.split(sbml)[1]
            if len(gbank.split(',')) == 1:
                gbank_dir, gbank = os.path.split(gbank)[0], os.path.split(gbank)[1]
            else:
                gbankout = []
                for g in gbank.split(','):
                    gbank_dir, gbanktmp = os.path.split(g)[0], os.path.split(g)[1]
                    gbankout.append(gbanktmp)
                gbank = gbankout
                del gbankout, gbanktmp

            modlist = [(key, gbank, sbml)]
            print(modlist)
            print(sbml_dir)
            print(gbank_dir)
            if sbml_dir != gbank_dir:
                QApplication.restoreOverrideCursor()
                title = 'Create SeQPlus Template Model: error'
                msg = 'Model creation error: input files must be in same directory'
                self.widgetMsgBox(QMessageBox.Critical, title, msg)
                return
            else:
                try:
                    ret = biotools.createSeqplusModel(modlist, sbml_dir, self.seqplus_files, 'usr',\
                                                      gene_db=self._genedb_path_,\
                                                      add_cobra_annot=False)
                except AssertionError:
                    QApplication.restoreOverrideCursor()
                    title = 'Create SeQPlus Template Model: error'
                    msg = 'Model creation error: Model must be encoded using SBML Level 3 with FBC version 1 or 2 format!'
                    self.widgetMsgBox(QMessageBox.Critical, title, msg)
                    return
                self.bp_lview.clear()
                datF = []
                for f in os.listdir(self.seqplus_files):
                    if f.endswith('.seqplus.xml'):
                        datF.append(f.replace('.seqplus.xml',''))
                datF.sort()
                self.bp_lview.addItems(datF)
                self.bp_lview.update()
        QApplication.restoreOverrideCursor()
        title = 'Create SeQPlus Template Model ({})'.format(ret[0])
        msg = 'Created new SeqPlus template model: ({})-({})'.format(ret[0], sbml)
        self.widgetMsgBox(QMessageBox.Information, title, msg)
        self.status_bar.showMessage(msg)
        self.widgetAddSeQPlus.setWindowTitle(title)
        self.widgetAddSeQPlus.activateWindow()

    @pyqtSlot()
    def menu_resetGUI(self):
        self.widget_tabpanel.setTabEnabled(0, True)
        self.widget_tabpanel.setCurrentIndex(0)
        try:
            for t in [3,2,1]:
                self.widget_tabpanel.removeTab(t)
        except Exception as e:
            print(e)
        self.menuSessions.setDisabled(True)
        #self.appwindow.setWindowTitle('MetaDraft')
        self.sbmlApp0.setDisabled(True)
        self.sbmlApp1.setDisabled(True)
        self.sbmlApp2.setDisabled(True)
        self.sbmlAppArch.setDisabled(True)
        self.excelApp.setDisabled(True)
        self.resetApp.setDisabled(True)
        self.tblxApp.setDisabled(True)
        self.menuTools.setEnabled(True)
        self.menuToolsM.setDisabled(True)
        self.widgetDisplayReact_update('Ready...')
        self.selected_metabolites = {}
        self.selected_reactions = {}
        self._gene_score_limits_ = [1.0, 1.0]
        #self.menuConfig.setEnabled(True)
        self._wnotes_.setDisabled(True)

    def func_getGroupMembership(self, mod):
        """
        Returns group membership of items in groups. Returns {object_id: ['group_id1', 'group_id2']}

        """
        grps = {}
        for g in mod.groups:
            gid = g.getId()
            for mid in g.getMemberIDs():
                if mid in grps:
                    grps[mid].append(gid)
                else:
                    grps[mid] = [gid]
        return grps

    def buildSBMLModel(self):
        self.menu_buildAll()
        self.model = cbmpy.CBModel.Model(self.model_name)
        user = self.func_getCurrentUser(realname=False)
        tmp_compartments = [self.default_compartment]
        for s in self.selected_metabolites:
            S = self.selected_metabolites[s].clone()
            #S.setAnnotation('metadraft_template', S._organism_)
            if S.getCompartmentId() not in tmp_compartments:
                tmp_compartments.append(S.getCompartmentId())
            #S.setCompartmentId(self.default_compartment)
            #S.unsetBoundary()
            notes = self.readNotesFromNotesDB(S._organism_, S.getId(), user, self._DAT_NOTESDB_KEY_)
            if notes != '':
                S.setNotes('<p>{}</p>'.format(notes))
            self.model.addSpecies(S)

        for c_ in tmp_compartments:
            self.model.createCompartment(c_)


        selected_reactions = ['{}{}{}'.format(str(self.table_reaction.item(ridx, 0).text()).strip(), self.id_sep,\
                                              str(self.table_reaction.item(ridx, 3).text()).strip())\
                              for ridx in range(self.table_reaction.rowCount()) if self.table_reaction.item(ridx, 2).checkState()]

        print('NSR', len(selected_reactions))
        print('OSR', len(self._reaction_selected_ids_))

        dup_re_db = {}

        gk = self._gene_selected_map_.keys()
        slx = []
        tlx = []
        new2old = {}
        for x in gk:
            x = x.split(self.id_sep)
            slx.append(x[0][len(self.gene_prefix):])
            tlx.append(x[1])
            new2old.update({x[1] : x[0][len(self.gene_prefix):]})
        #print(slx)
        #print(tlx)
        #print(new2old)

        new_groups = {}
        for r_ in selected_reactions:
            R = self.selected_reactions[r_]['obj'].clone()
            #R.setAnnotation('metadraft_template', R._organism_)

            #print('-', R.__objref__)
            #for r in R.reagents:
                #print('--', r.__objref__)

            rid = R.getId()

            # store group information
            if R._organism_ not in new_groups:
                new_groups[R._organism_] = {'grp_map' : self.func_getGroupMembership(self._DAT_MODELS[R._organism_]),
                                            'grpd' : [rid]}
            else:
                new_groups[R._organism_]['grpd'].append(rid)


            R.setCompartmentId(self.default_compartment)
            DUP_R = False
            notes = self.readNotesFromNotesDB(R._organism_, R.getId(), user, self._DAT_NOTESDB_KEY_)
            if notes != '':
                R.setNotes('<p>{}</p>'.format(notes))
            try:
                self.model.addReaction(R, create_default_bounds=False)
                self.model.createReactionBounds(rid, self.selected_reactions[r_]['obj'].getLowerBound(), self.selected_reactions[r_]['obj'].getUpperBound())
                ridnew = rid
            except RuntimeError:
                #R.__objref__ = None
                if rid in dup_re_db:
                    dup_re_db[rid] = dup_re_db[rid] + 1
                else:
                    dup_re_db[rid] = 1
                ridnew = '{}_copy_{}'.format(rid, dup_re_db[rid])
                dup_re_db[ridnew] = dup_re_db[rid]
                #print(dup_re_db)
                R.setId(ridnew)
                self.model.addReaction(R, create_default_bounds=False)
                # add equation as
                self.model.createReactionBounds(ridnew, self.selected_reactions[r_]['obj'].getLowerBound(), self.selected_reactions[r_]['obj'].getUpperBound())
                print('Suspected duplicate reaction ID detected {} included as {}!'.format(rid, ridnew))

            gpr = self._DAT_MODELS[R._organism_].getGPRforReaction(rid)
            gene_ass = gene_ass_original = gpr.getAssociationStr()

            gids = gpr.getGeneLabels()
            gids.sort()
            gids.reverse()

            used_gids = []
            altlabels = {}
            #print(self._gene_selected_ids_)
            for g_ in gids:
                G = self._DAT_MODELS[R._organism_].getGeneByLabel(g_)
                if g_ not in self._gene_selected_ids_:
                    #gene_ass = gene_ass.replace(g_, 'unknown')
                    altlabels[G.getId()] = 'unknown'
                else:
                    used_gids.append(g_)
                    altlabels[G.getId()] = new2old[G.getLabel()]

            #for g_ in used_gids:
                #altlabels[G.getId()] = slx[tlx.index(g_)]

            if not DUP_R:
                self.model.createGeneProteinAssociation(ridnew, gene_ass, altlabels=altlabels)
                R.setAnnotation('GENE_ASSOCIATION_ORIGINAL', gene_ass_original)
                if self.NO_EXPORT_SEQ:
                    for a_ in list(R.annotation):
                        if a_.startswith('gbank_seq_'):
                            R.annotation.pop(a_)

        del_unknowns = []
        for g in self.model.genes:
            if g.getLabel() == 'unknown':
                del_unknowns.append(g.getId())

        _gene_selected_labels_ = [self.model.getGeneIdFromLabel(new2old[g]) for g in self._gene_selected_ids_]
        for GPR in self.model.gpr:
            for gid in GPR.getGeneIds():
                if gid not in _gene_selected_labels_:
                    if gid not in del_unknowns:
                        del_unknowns.append(gid)
                        print('Added nasty gene:', gid)


        for gid in list(set(del_unknowns)):
            self.model.deleteGene(gid)

        #self.model.deleteGene('G_b2416')

        for g in self.model.genes:
            g.setId('{}_{}'.format(cbmpy.CBModel.fixId(g.getLabel()), g.getId()))

        for gpr in self.model.gpr:
            self.model.getReaction(gpr.getProtein()).setAnnotation('GENE_ASSOCIATION', gpr.getAssociationStr())

        model_gene_ids = self.model.getGeneIds()
        for g_ in model_gene_ids:
            notes = self.readNotesFromNotesDB(g_, g_, user, self._DAT_NOTESDB_KEY_)
            if notes != '':
                self.model.getGene(g_).setNotes('<p>{}</p>'.format(notes))
                #print('gn', self.model.getGene(g_).getNotes())

        model_gene_lbls = self.model.getGeneLabels()
        #for r in range(self.table_gene.rowCount()):
            #srcg = str(self.table_gene.item(r, 0).text())[len(self.gene_prefix):]
            #matchg = str(self.table_gene.item(r, 1).text())
            #score = str(self.table_gene.item(r, 2).text())
            #org = str(self.table_gene.item(r, 4).text())
            ##print(srcg, matchg, score, org)
            #if srcg in model_gene_lbls:
                #G = self.model.getGeneByLabel(srcg)
                #G.setAnnotation('metadraft_match', matchg)
                #G.setAnnotation('metadraft_score', score)
                #G.setAnnotation('metadraft_template', org)

        # generate new group information
        cntr = 1
        grps = {}

        for o in new_groups:
            grps = {}
            tmpmod = self._DAT_MODELS[o]
            for mid in new_groups[o]['grpd']:
                if mid in new_groups[o]['grp_map']:
                    for gid in new_groups[o]['grp_map'][mid]:
                        if gid not in grps:
                            grps[gid] = []
                        grps[gid].append(mid)

            for g in grps:
                new_gid = '{}_{}'.format(cbmpy.CBModel.fixId(o), g)
                self.model.createGroup(new_gid)
                G = self.model.getGroup(new_gid)
                G.setName(tmpmod.getGroup(g).getName())
                for mid in grps[g]:
                    # a bit of a hack for now, only deals with reactions and species
                    obj = self.model.getReaction(mid)
                    if obj is None:
                        obj = self.model.getSpecies(mid)
                    if obj is not None:
                        G.addMember(obj)

            new_gid = 'metadraft_{}'.format(cbmpy.CBModel.fixId(o))
            self.model.createGroup(new_gid)
            G = self.model.getGroup(new_gid)
            G.setName('MetaDraft reactions {}'.format(o))
            G.setNotes('These {} reactions were included from the {} model using MetaDraft (https://github.com/SystemsBioinformatics/metadraft) ver. {}'.format(len(new_groups[o]['grpd']), o, __version__))
            G.addMember([self._DAT_MODELS[o].getReaction(r_) for r_ in new_groups[o]['grpd']])

        return self.model

    @pyqtSlot()
    def menu_exportSBML0(self):
        self.menu_exportSBML(fbcv=0)

    @pyqtSlot()
    def menu_exportSBML1(self):
        self.menu_exportSBML(fbcv=1)

    @pyqtSlot()
    def menu_exportSBML2(self):
        self.menu_exportSBML(fbcv=2)

    @pyqtSlot()
    def menu_exportCOMBINE(self):

        arch_fname = str(self.saveFile('Export COMBINE archive', self.sbml_save_dir, '*.omex'))
        print(arch_fname)
        out_dir, arch_fname = os.path.split(arch_fname)
        fname = arch_fname.replace('.omex', '')
        print(out_dir, fname)
        #fname = 'combine-export'
        print('File: {}'.format(fname))
        print('Dir: {}'.format(self._tmpDir_))


        self.widget_busy.setWindowTitle("Exporting COMBINE archive")
        self.widget_busy.setLabelText("Exporting COMBINE archive")
        self.widget_busy.show()
        self.widgetBusyUpdate(10)

        arch_dir = os.path.join(self._tmpDir_, '_combine_')
        if os.path.exists(arch_dir):
            for f_ in os.listdir(arch_dir):
                print(f_)
                os.remove(os.path.join(arch_dir, f_))
        else:
            os.makedirs(arch_dir)

        model = self.buildSBMLModel()
        if not fname.endswith('.xml'):
            fname += '.xml'

        for r_ in model.reactions:
            eqn = r_.getEquation()
            if eqn is not None:
                r_.setAnnotation('equation_id', eqn)
            eqn = r_.getEquation(use_names=True)
            if eqn is not None:
                r_.setAnnotation('equation_name', eqn)

        zfpath = os.path.join(out_dir, fname.replace('.xml', '.omex'))
        print(zfpath)
        zf = zipfile.ZipFile(zfpath, mode='w', compression=zipfile.ZIP_DEFLATED)

        self.widgetBusyUpdate(20)
        cbmpy.writeSBML3FBCV2(model, os.path.join(arch_dir, fname), add_cobra_annot=False)
        zf.write(os.path.join(arch_dir, fname), arcname=fname)
        self.widgetBusyUpdate(30)
        cbmpy.writeModelToExcel97(model, os.path.join(arch_dir, fname.replace('.xml', '')))
        zf.write(os.path.join(arch_dir, fname.replace('.xml', '.xls')), arcname=fname.replace('.xml', '.xls'))
        self.widgetBusyUpdate(40)
        self.func_generateSummaryReport()
        self.widgetBusyUpdate(50)
        srpt = self.func_formatSummaryReport()
        self.widgetBusyUpdate(60)
        grpt = self.func_formatGeneReport()
        self.widgetBusyUpdate(70)
        rrpt = self.func_formatReactionReport()
        self.widgetBusyUpdate(80)
        mrpt = self.func_formatMetaboliteReport()
        self.widgetBusyUpdate(90)

        tfn = os.path.join(arch_dir, '1_summary_report.html')
        F = open(tfn, 'w')
        F.write(srpt)
        F.close()
        zf.write(tfn, arcname='1_summary_report.html')

        tfn = os.path.join(arch_dir, '2_gene_report.html')
        F = open(tfn, 'w')
        F.write(grpt)
        F.close()
        zf.write(tfn, arcname='2_gene_report.html')

        tfn = os.path.join(arch_dir, '3_reaction_report.html')
        F = open(tfn, 'w')
        F.write(rrpt)
        F.close()
        zf.write(tfn, arcname='3_reaction_report.html')

        tfn = os.path.join(arch_dir, '4_metabolite_report.html')
        F = open(tfn, 'w')
        F.write(mrpt)
        F.close()
        zf.write(tfn, arcname='4_metabolite_report.html')

        self.widgetBusyUpdate(95)
        F = open(os.path.join(arch_dir, 'index.html'), 'w')
        F.write(metadraftreports.combine_index_template(fname, fname.replace('.xml', '.xls')))
        F.close()
        zf.write(os.path.join(arch_dir, 'index.html'), arcname='index.html')

        F = open(os.path.join(arch_dir, 'manifest.xml'), 'w')
        F.write(metadraftreports.combine_manifest_file(fname, fname.replace('.xml', '.xls')))
        F.close()
        zf.write(os.path.join(arch_dir, 'manifest.xml'), arcname='manifest.xml')

        F = open(os.path.join(arch_dir, 'metadata.rdf'), 'w')
        F.write(metadraftreports.combine_metadata_file(vc_given=self.func_getCurrentUser(),
                                                      vc_family='MetaDraft Software',
                                                      vc_email='', vc_org='', scTime=None))
        F.close()
        zf.write(os.path.join(arch_dir, 'metadata.rdf'), arcname='metadata.rdf')

        zf.close()
        self.widgetBusyUpdate(100)


    @pyqtSlot()
    def menu_exportSBML(self, fbcv):
        if fbcv == 0:
            filename = str(self.saveFile('Export COBRA SBML dialect', self.sbml_save_dir, '*.xml'))
        else:
            filename = str(self.saveFile('Export SBML L3 FBCv{}'.format(fbcv), self.sbml_save_dir, '*.xml'))
        self.widget_busy.setWindowTitle("Exporting SBML")
        self.widget_busy.setLabelText("Exporting SBML")
        self.widget_busy.show()
        self.widgetBusyUpdate(10)
        if filename == '' or filename == None:
            self.widget_busy.cancel()
            return
        elif not filename.endswith('.xml'):
            filename += '.xml'
        model = self.buildSBMLModel()
        self.widgetBusyUpdate(50)

        fsplit = os.path.split(filename)
        self.sbml_save_dir = fsplit[0]
        self.sbml_file_name = fsplit[1]
        if fbcv == 0:
            cbmpy.writeCOBRASBML(model, filename)
        elif fbcv == 1:
            cbmpy.writeSBML3FBC(model, filename, add_cobra_annot=True)
        elif fbcv == 2:
            cbmpy.writeSBML3FBCV2(model, filename, add_cobra_annot=False)
        self.sbml_model = model
        self.widgetBusyUpdate(100)
        del model
        self.status_bar.showMessage('SBML export successful: {}'.format(filename))

    @pyqtSlot()
    def menu_exportExcel(self):
        filename = str(self.saveFile('Export Excel spreadsheet', self.excel_save_dir, '*.xls'))
        self.widget_busy.setWindowTitle("Exporting Model to Excel")
        self.widget_busy.setLabelText("Exporting Model to Excel")
        self.widget_busy.show()
        self.widgetBusyUpdate(10)
        if filename == '' or filename == None:
            self.widget_busy.cancel()
            return
        elif not filename.endswith('.xls'):
            filename += '.xls'

        model = self.buildSBMLModel()

        self.widgetBusyUpdate(40)

        for r_ in model.reactions:
            eqn = r_.getEquation()
            if eqn is not None:
                r_.setAnnotation('equation_id', eqn)
            eqn = r_.getEquation(use_names=True)
            if eqn is not None:
                r_.setAnnotation('equation_name', eqn)

        self.widgetBusyUpdate(50)

        fsplit = os.path.split(filename)
        self.excel_save_dir = fsplit[0]
        self.excel_file_name = fsplit[1]
        cbmpy.writeModelToExcel97(model, filename[:-4])
        self.sbml_model = model
        self.widgetBusyUpdate(100)
        del model
        self.status_bar.showMessage('Excel export successful: {}'.format(filename))

    @pyqtSlot()
    def menu_exportTableToCSV(self):
        if self._tabpanel_idx_[self._active_tab_] == 'Genes':
            self.func_tableSaveCSV(self.table_gene)
        elif self._tabpanel_idx_[self._active_tab_] == 'Reactions':
            self.func_tableSaveCSV(self.table_reaction)
        elif self._tabpanel_idx_[self._active_tab_] == 'Metabolites':
            self.func_tableSaveCSV(self.table_metab)

    @pyqtSlot()
    def menu_convertSBML(self):
        print("menu_convertSBML")

        self.widget_menu_sbmlConvertApp = QWidget()
        self.widget_menu_sbmlConvertApp.setWindowTitle('SBML tools')
        ppos = self.func_getNewPopupWindowCoords()
        self.widget_menu_sbmlConvertApp.setGeometry(QtCore.QRect(ppos[0], ppos[1], 450, 150))


        fnamel = QLabel(self.widget_menu_sbmlConvertApp)
        fnamel.setText('File name')
        ftypel = QLabel(self.widget_menu_sbmlConvertApp)
        ftypel.setText('File type')
        fname = QLineEdit(self.widget_menu_sbmlConvertApp)
        fname.setReadOnly(True)
        ftype = QLineEdit(self.widget_menu_sbmlConvertApp)
        ftype.setReadOnly(True)

        @pyqtSlot()
        def openFile():
            F = None
            output = None
            try:
                F = os.path.abspath(str(self.openFile('SBML file', self._history_open_dir_, '*')))
                fname.setText(F)
                output = msg = ''
                try:
                    if F is not None:
                        output, msg = cbmpy.CBXML.sbml_fileFindVersion(F)
                except IOError:
                    print('ERROR: SBML convert: invalid filename')
                ftype.setText('{} ({})'.format(output, msg.split(',')[0]))
            except IOError:
                print('ERROR: SBML convert: invalid filename')
            if output == 'L3V1FBC1':
                toFBC1.setEnabled(False)
                toFBC2.setEnabled(True)
            elif output == 'L3V1FBC2':
                toFBC1.setEnabled(True)
                toFBC2.setEnabled(False)
            elif output == 'L2FBA' or output == 'COBRA':
                toFBC1.setEnabled(True)
                toFBC2.setEnabled(True)
            self._test_sbml_file_ = F
            self._test_sbml_type_ = output
            self.widget_menu_sbmlConvertApp.activateWindow()

        openbut = QPushButton(self.widget_menu_sbmlConvertApp)
        openbut.setText('Load SBML')
        openbut.clicked.connect(openFile)

        #typebut = QPushButton(self.widget_menu_sbmlConvertApp)
        #typebut.setText('Get Type')
        #typebut.clicked.connect(getType)

        @pyqtSlot()
        def toFBC1F():
            xmod = None
            #try:
            if self._test_sbml_type_ == 'COBRA':
                xmod = cbmpy.readCOBRASBML(self._test_sbml_file_)
            elif self._test_sbml_type_ == 'L2FBA':
                xmod = cbmpy.readSBML2FBA(self._test_sbml_file_)
            elif self._test_sbml_type_ == 'L3V1FBC1' or self._test_sbml_type_ == 'L3V1FBC2':
                xmod = cbmpy.readSBML3FBC(self._test_sbml_file_)
            if xmod is None:
                raise RuntimeError('Invalid file to convert')
            #except:
                #xmod = None
                #self.widgetMsgBox(QMessageBox.Critical, 'SBML Conversion Error', 'File conversion failed!')
            if xmod is not None:
                cbmpy.writeSBML3FBC(xmod, str(self.saveFile('Save converted file', self._history_save_dir_, '*.fbcv1.xml')))
            del xmod

        @pyqtSlot()
        def toFBC2F():
            xmod = None
            #try:
            if self._test_sbml_type_ == 'COBRA':
                xmod = cbmpy.readCOBRASBML(self._test_sbml_file_)
            elif self._test_sbml_type_ == 'L2FBA':
                xmod = cbmpy.readSBML2FBA(self._test_sbml_file_)
            elif self._test_sbml_type_ == 'L3V1FBC1' or self._test_sbml_type_ == 'L3V1FBC2':
                xmod = cbmpy.readSBML3FBC(self._test_sbml_file_)
            if xmod is None:
                raise RuntimeError('Invalid file to convert')
            #except:
                #xmod = None
                #self.widgetMsgBox(QMessageBox.Critical, 'SBML Conversion Error', 'File conversion failed!')
            if xmod is not None:
                cbmpy.writeSBML3FBCV2(xmod, str(self.saveFile('Save converted file', self._history_save_dir_, '*.fbcv2.xml')))
            del xmod


        toFBC1 = QPushButton(self.widget_menu_sbmlConvertApp)
        toFBC1.setText('Save as SBML3 FBCv1')
        toFBC1.setEnabled(False)
        toFBC1.clicked.connect(toFBC1F)

        toFBC2 = QPushButton(self.widget_menu_sbmlConvertApp)
        toFBC2.setText('Save as SBML3 FBCv2')
        toFBC2.setEnabled(False)
        toFBC2.clicked.connect(toFBC2F)


        layout = QGridLayout(self.widget_menu_sbmlConvertApp)
        layout.setSpacing(10)
        layout.addWidget(fnamel, 0, 0, 1, 2)
        layout.addWidget(fname, 0, 1, 1, 3)
        layout.addWidget(ftypel, 1, 0, 1, 3)
        layout.addWidget(ftype, 1, 1, 1, 3)
        layout.addWidget(openbut, 2, 0)
        #layout.addWidget(typebut, 2, 1)
        layout.addWidget(toFBC1, 2, 1)
        layout.addWidget(toFBC2, 2, 2)

        self.widget_menu_sbmlConvertApp.show()

    def func_tableSaveCSV(self, table):
        path = self.saveFile('Save File', directory=self._history_save_dir_, filterex='CSV(*.csv)')
        if path is not None and len(path) > 0:
            if not str(path).endswith('.csv'):
                path += '.csv'
            with open(unicode(path), 'wb') as stream:
                writer = csv.writer(stream)
                for row in range(table.rowCount()):
                    rowdata = []
                    for column in range(table.columnCount()):
                        item = table.item(row, column)
                        if item.checkState() == Qt.Checked:
                            rowdata.append(unicode('checked').encode('utf8'))
                        # this can be removed when the gene prefix is killed
                        elif item is not None and column == 0 and self._tabpanel_idx_[self._active_tab_] == 'Genes':
                            rowdata.append(unicode(item.text()).encode('utf8')[len(self.gene_prefix):])
                        elif item is not None:
                            rowdata.append(unicode(item.text()).encode('utf8'))
                        else:
                            rowdata.append('')
                    writer.writerow(rowdata)

    def func_findNonGprReactionsForModel(self, model):
        # get all reaction id's
        all_reac_ids = model.getReactionIds()

        # get gene associated reaction id's
        gene_assoc_reac_ids = []
        for gpr in model.gpr:
            gene_assoc_reac_ids.append(gpr.getProtein())

        # the difference between these two sets are the non-gene associated reactions
        non_gpr_reac_ids = list(set(all_reac_ids).difference(set(gene_assoc_reac_ids)))

        ## now we want to write out a csv with the ids, names and equations of the non-gpr reactions
        #for rid in non_gpr_reac_ids:
            #R = model.getReaction(rid)
            #rows.append([rid, R.getName(), R.getEquation(), R.getEquation(use_names=True)[:-1]])

        return non_gpr_reac_ids

    def widgetBusy(self):
        self.widget_busy = QProgressDialog(self)
        #self.widget_busy.setWindowModality(Qt.WindowModal)
        self.widget_busy.setRange(0,100)
        self.widget_busy.setMinimumDuration(0)
        self.widget_busy.setValue(100)

    def widgetBusyUpdate(self, value):
        self.widget_busy.setValue(value-1)
        QCoreApplication.instance().processEvents()
        self.widget_busy.setValue(value)
        QCoreApplication.instance().processEvents()

    def saveFile(self, title, directory, filterex):
        filename = QFileDialog.getSaveFileName(self, title, directory, filter=filterex)
        if HAVE_QT5:
            filename = filename[0]
        self._history_save_dir_ = os.path.dirname(str(filename))
        return str(filename)

    def openFile(self, title, directory, filterex):
        filename = QFileDialog.getOpenFileName(self, title, directory, filter=filterex)
        if HAVE_QT5:
            filename = filename[0]
        filename = os.path.abspath(str(filename))
        self._history_open_dir_ = os.path.dirname(filename)
        return str(filename)

    def openFiles(self, title, directory, filterex):
        filenames = QFileDialog.getOpenFileNames(self, title, directory, filter=filterex)
        if HAVE_QT5:
            filenames = [os.path.abspath(str(a[0])) for a in filenames]
        else:
            filenames = [os.path.abspath(str(a)) for a in filenames]
        self._history_open_dir_ = os.path.dirname(str(filenames[0]))
        return filenames

    def widgetDisplayReact(self):
        self.reactDisplay = QTextBrowser()
        self.reactDisplay.setOpenExternalLinks(True)
        self.reactDisplay.setHtml(str('<html>Ready</html>'))

    def widgetDisplayReact_dev(self):
        self.reactDisplay = QtWebKit.QWebView(self)
        self.reactDisplay.setHtml(str('<html>Ready</html>'))

    def widgetDisplayReact_update(self, txt):
        self.reactDisplay.setHtml(str(txt))
        self.reactDisplay.update()

    def widgetTableReaction(self, parent=None):
        self.table_reaction = QTableWidget(0, len(self.table_react_cols), parent=self.widget_tabpanel)
        self.table_reaction.setColumnWidth(0,100)
        self.widgetTableReaction_populate()
        self.table_reaction.setSortingEnabled(True)

        # this works for mouseclicks
        #self.table_reaction.cellClicked.connect(self.widgetTableReaction_cellClicked)
        # this should also work for cursor changes
        self.table_reaction.itemSelectionChanged.connect(self.widgetTableReaction_cellSelectionChanged)
        self.table_reaction.cellChanged.connect(self.widgetTableReaction_cellChecked)
        # not sure if I want this right now until I have sorted out the persistence between sessions issues
        #self.table_reaction.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.table_reaction.customContextMenuRequested.connect(self.widgetTableReaction_tableRightClicked)

    def widgetTableReaction_populate(self):
        self.reaction_table_loading = True
        self.table_reaction.clear()
        self.table_reaction.setRowCount(0)
        self.table_reaction.setColumnCount(0)
        selected_reactions = {}
        for row in range(self.table_gene.rowCount()):
            if self.table_gene.item(row, 3).checkState() == QtCore.Qt.Checked:
                geneId = str(self.table_gene.item(row, 1).text())
                new_geneId = str(self.table_gene.item(row, 0).text())[len(self.gene_prefix):]
                geneIdtarget = str(self.table_gene.item(row, 0).text())
                if geneId in self._DAT_G2REACT:
                    for r_ in self._DAT_G2REACT[geneId]:
                        rid = str(r_.getId())
                        org = str(self._DAT_LINK_DICT_['__idx__'][geneId])
                        rkey = '{}{}{}'.format(rid, self.id_sep, org)
                        r_._organism_ = org
                        if rkey in selected_reactions:
                            selected_reactions[rkey]['genes'].append(geneIdtarget)
                            selected_reactions[rkey]['genematch'].append(geneId)
                            selected_reactions[rkey]['genematch_new'].append(new_geneId)
                        else:
                            selected_reactions[rkey] = {'id' : rid,
                                                       'name' : r_.getName(),
                                                       'org' : org,
                                                       'genes' : [geneIdtarget],
                                                       'genematch' : [geneId],
                                                       'genematch_new' : [new_geneId],
                                                       'obj' : r_
                                                       }
                else:
                    print('WARNING: {} not in _DAT_G2REACT: gene defined that does not appear in a GPR association').format(geneId)
            else:
                geneId = str(self.table_gene.item(row, 1).text())
                if geneId in self._DAT_G2REACT:
                    for r_ in self._DAT_G2REACT[geneId]:
                        r_._organism_ = str(self._DAT_LINK_DICT_['__idx__'][geneId])
        self.table_reaction.setRowCount(len(selected_reactions))
        self.table_reaction.setColumnCount(len(self.table_react_cols))
        self.table_reaction.setHorizontalHeaderLabels(self.table_react_cols)
        row = 0
        sel_reac = list(selected_reactions.keys())
        sel_reac.sort
        for r_ in sel_reac:
            item0 = QTableWidgetItem(selected_reactions[r_]['id'])
            item1 = QTableWidgetItem(selected_reactions[r_]['name'])
            item2 = QTableWidgetItem()
            item3 = QTableWidgetItem(selected_reactions[r_]['org'])
            # optional, show genematch the alternative is the genesource
            #item4 = QTableWidgetItem(','.join(selected_reactions[r_]['genes']))
            item4 = NumberTableListLengthItem(','.join(selected_reactions[r_]['genematch']))
            item5 = NumberTableListLengthItem(','.join(selected_reactions[r_]['genematch_new']))

            items = [item0, item1, item2, item3, item4, item5]
            for c_ in range(len(items)):
                if c_ == 2:
                    items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                    items[c_].setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                    if r_ in self._reaction_selected_map_:
                        #print('INFO: duplicate reaction ID in reaction table')
                        if self._reaction_selected_map_[r_]:
                            items[c_].setCheckState(QtCore.Qt.Checked)
                        else:
                            items[c_].setCheckState(QtCore.Qt.Unchecked)
                    else:
                        items[c_].setCheckState(QtCore.Qt.Checked)
                else:
                    items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                self.table_reaction.setItem(row, c_, items[c_])
            row += 1
        self.selected_reactions = selected_reactions
        #if len(self._reaction_selected_map_) > 0:
        #self._reaction_selected_map_ = self.widgetTableReaction_getMap()
        #self._reaction_selected_ids_ = self.widgetTableReaction_getSelectedIds()
        self.reaction_table_loading = False
        self._updateReactionMap_()
        self.table_reaction.sortByColumn(0, Qt.SortOrder(0))


    @pyqtSlot(QtCore.QPoint)
    def widgetTableReaction_tableRightClicked(self, QPos):
        print('REACTION TABLE RIGHT CLICKED', QPos)

        self.widget_tablereact_rclickmenu = QMenu()
        menu_item = self.widget_tablereact_rclickmenu.addAction("Edit annotation")
        menu_item.triggered.connect(self.menu_react_editAnnotation)
        parentPosition = self.table_reaction.viewport().mapToGlobal(QPos)
        self.widget_tablereact_rclickmenu.move(parentPosition)
        self.widget_tablereact_rclickmenu.show()

    @pyqtSlot()
    def menu_react_editAnnotation(self):
        idx = self.table_reaction.selectedIndexes()[0]
        self.table_reaction.selectRow(idx.row())
        item = self.table_reaction.item(idx.row(), 0)
        rid = item.text()
        mod = self._DAT_MODELS[str(self.table_reaction.item(idx.row(), 3).text())]
        self.widget_annot_edit = CBMPyAnnotationEditor(mod, rid)
        parentPosition = self.table_reaction.viewport().mapToGlobal(self.pos())
        parentPosition.setY(parentPosition.y()-250)
        self.widget_annot_edit.move(parentPosition)

    def widgetBuildPanel(self):
        bp = QWidget(parent=self.widget_tabpanel)
        bp.layout = QGridLayout()
        bp.layout.setSpacing(10)

        bp_label_targ = QLabel(bp)
        bp_label_targ.setText('User proteome (sequence)')
        bp_btn_targ = QPushButton(bp)
        bp_btn_targ.setText('Select')
        bp_btn_targ.clicked.connect(self.bp_targetOpen)
        bp_text_targ = QLineEdit(bp)
        bp_text_targ.setReadOnly(True)

        bp_label_out = QLabel(bp)
        bp_label_out.setText('Benchmark (FASTA)')
        bp_label_out.setDisabled(True)
        bp_btn_out = QPushButton(bp)
        bp_btn_out.setText('Select')
        bp_btn_out.setDisabled(True)
        bp_btn_out.clicked.connect(self.bp_targetOpen2)
        bp_text_out = QLineEdit(bp)
        bp_text_out.setDisabled(True)

        bp_label_lib = QLabel(bp)
        bp_label_lib.setText('Select target network(s) ')

        bp_btn_meta = QPushButton(bp)
        bp_btn_meta.setText('Build metaproteome')
        bp_btn_meta.setEnabled(False)
        bp_btn_meta.clicked.connect(self.bp_buildMetaproteome)

        bp_btn_blast = QPushButton(bp)
        bp_btn_blast.setText('BLAST')
        bp_btn_blast.setEnabled(False)
        bp_btn_blast.clicked.connect(self.bp_runBLAST)

        bp_lview = QListWidget(bp)
        bp_lview.setDragDropMode(QAbstractItemView.InternalMove)
        bp_lview.setSelectionMode(QAbstractItemView.ExtendedSelection)
        datF = []
        for f in os.listdir(self.seqplus_files):
            if f.endswith('.seqplus.xml'):
                datF.append(f.replace('.seqplus.xml',''))
        datF.sort()
        bp_lview.addItems(datF)
        bp_lview.setContextMenuPolicy(Qt.CustomContextMenu)
        bp_lview.customContextMenuRequested.connect(self.bp_lviewRightClicked)
        bp_label_res = QLabel(bp)
        bp_label_res.setText('Load result')

        ## this is the simple directory listing
        #bp_lview_res = QListWidget(bp)
        #bp_lview_res.setDragDropMode(QAbstractItemView.InternalMove)
        #bp_lview_res.setSelectionMode(QAbstractItemView.ExtendedSelection)
        #bp_lview_res.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #bp_lview_res.customContextMenuRequested.connect(self.bp_lviewOnRightClick)
        #datF = []
        #for f in os.listdir(self.result_files):
            #if f.endswith('.resplus.json'):
                #datF.append(f.replace('.resplus.json',''))
        #datF.sort()
        #bp_lview_res.addItems(datF)

        ## this is the replacement tree view for the result listing
        self.widgetResultTree()

        bp_btn_res = QPushButton(bp)
        bp_btn_res.setText('Load')
        bp_btn_res.clicked.connect(self.bp_loadResult)

        bp.layout.addWidget(bp_label_targ, 1, 0)
        bp.layout.addWidget(bp_text_targ, 1, 1)
        bp.layout.addWidget(bp_btn_targ, 1, 2)

        bp.layout.addWidget(bp_label_out, 2, 0)
        bp.layout.addWidget(bp_text_out, 2, 1)
        bp.layout.addWidget(bp_btn_out, 2, 2)

        bp.layout.addWidget(bp_label_lib, 3, 0)
        bp.layout.addWidget(bp_lview, 3, 1, 4, 1)
        bp.layout.addWidget(bp_btn_meta, 3, 2)
        bp.layout.addWidget(bp_btn_blast, 5, 2)
        bp.layout.addWidget(bp_label_res, 8, 0)

        ## simple view
        #bp.layout.addWidget(bp_lview_res, 6, 1, 2, 1)

        # tree view
        bp.layout.addWidget(self.tree_results, 8, 1, 2, 1)
        bp.layout.addWidget(bp_btn_res, 8, 2)
        # replaced with menu option
        #bp.layout.addWidget(self.bp_btn_idopt, 4, 0)
        bp.setLayout(bp.layout)

        self.build_panel = bp
        self.bp_text_targ = bp_text_targ
        self.bp_text_out = bp_text_out
        self.bp_label_out = bp_label_out
        self.bp_btn_out = bp_btn_out
        self.bp_lview = bp_lview
        #self.bp_lview_res = bp_lview_res
        self.bp_btn_meta = bp_btn_meta
        self.bp_btn_blast = bp_btn_blast

    def bp_targetOpen(self):
        self.widget_busy.setWindowTitle("Loading sequences")
        self.widget_busy.setLabelText("Loading sequences")
        self.widget_busy.show()
        fname = self.openFile('Open Target file', self._history_open_dir_, 'Supported (*.gbk *.gbff *.gb *.fasta *.faa *.fa);;FASTA (*.fasta *.faa *.fa);;GenBank (*.gbk *.gbff *.gb)')
        self.widgetBusyUpdate(10)
        if fname.endswith('.fasta') or fname.endswith('.faa') or fname.endswith('.fa') or fname.endswith('.gbk') or fname.endswith('.gb') or fname.endswith('.gbff'):
            fname1 = biotools.createParanoidFASTAfromFile(fname, ext_replace='.in.fasta', gene_prefix=self.gene_prefix)
            self.widgetBusyUpdate(50)
            if not (fname.endswith('.fasta') or fname.endswith('.faa') or fname.endswith('.fa')):
                biotools.addGeneInformationToDB(fname, self._genedb_, 'GENES', 'id')
            self.widgetBusyUpdate(80)
            fname = fname1
            del fname1
        else:
            self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fname)))
            self.widgetBusyUpdate(100)
            return
        self.bp_text_targ.setText(fname)
        r_html = self.buildHtmlStringMetaprot()
        self.widgetBusyUpdate(90)
        self.widgetDisplayReact_update(r_html)
        self.status_bar.showMessage('User proteome: \"{}\" opened.'.format(str(self.bp_text_targ.text())))
        self.widgetBusyUpdate(100)
        self.bp_btn_meta.setEnabled(True)

    def bp_targetOpen2(self):
        fname = self.openFile('Open benchmark file', self._history_open_dir_, 'Supported (*.gbk *.gbff *.gb *.fasta *.faa *.fa);;FASTA (*.fasta *.faa *.fa);;GenBank (*.gbk *.gbff *.gb)')
        if fname.endswith('.fasta') or fname.endswith('.faa') or fname.endswith('.fa') or fname.endswith('.gbk') or fname.endswith('.gb') or fname.endswith('.gbff'):
            fname = biotools.createParanoidFASTAfromFile(fname)
        else:
            self.status_bar.showMessage('Invalid file type: \"{}\" ignored.'.format(str(fname)))
            return
        self.bp_text_out.setText(fname)
        r_html = self.buildHtmlStringMetaprot()
        self.widgetDisplayReact_update(r_html)
        self.status_bar.showMessage('Benchmark proteome: \"{}\" loaded.'.format(str(self.bp_text_out.text())))

    def bp_loadResult(self):
        #if len(self.bp_lview_res.selectedItems()) < 1:
            #return
        #print(self.result_file)
        if not os.path.exists(self.result_file):
            self.status_bar.showMessage('Invalid selection: \"{}\".'.format(self.result_file))
            return

        self.widget_busy.setWindowTitle("Loading results")
        self.widget_busy.setLabelText("Loading results")
        self.widget_busy.show()
        self.widgetBusyUpdate(10)
        self._colourCycler_ = itertools.cycle(self._colours_)
        self.widgetBusyUpdate(50)
        self.getMetaProteomeData(self.result_file)
        self.widgetBusyUpdate(70)
        self.widgetPanel()
        self.menuSessions.setEnabled(True)
        self.menu_savedSessionMenu()
        self.widgetBusyUpdate(80)
        self.widget_tabpanel.setTabEnabled(1, True)
        self.widget_tabpanel.setTabEnabled(2, True)
        self.widget_tabpanel.setTabEnabled(3, True)
        self.widget_tabpanel.setCurrentIndex(1)
        self.widget_tabpanel.setTabEnabled(0, False)
        self.notes = str(self._wnotes_.toHtml())
        self.widgetBusyUpdate(90)
        self.widgetResultTree_fill()
        self.sbmlApp0.setEnabled(True)
        self.sbmlApp1.setEnabled(True)
        self.sbmlApp2.setEnabled(True)
        self.sbmlAppArch.setEnabled(True)
        self.excelApp.setEnabled(True)
        self.resetApp.setEnabled(True)
        self.tblxApp.setEnabled(True)
        self.menuTools.setDisabled(True)
        self.menuToolsM.setEnabled(True)
        #self.menuConfig.setDisabled(True)
        self.widgetBusyUpdate(100)
        self.status_bar.showMessage('Loaded results: \"{}\".'.format(self.result_file))
        self.appwindow.setWindowTitle('MetaDraft {} result: {}'.format(__version__, os.path.split(self.result_file)[-1].replace('_metalink.resplus.json','')))


    def bp_runBLAST(self):
        if str(self.bp_text_targ.text()) == '' and self.metaproteome_file == '':
            self.status_bar.showMessage('ERROR: please select both a user target file and generate a metaproteome')
            self.widgetDisplayReact_update('<html><body><h3>ERROR: please input a user proteome and generate a metaproteome</h3></body></html>')
            return
        self.runBLASTFUNC()
        self.bp_btn_blast.setEnabled(False)
        self.bp_btn_blast.update()
        self.status_bar.showMessage('BLAST run successful')
        self.widgetResultTree_fill()

        #datF = []
        #for f in os.listdir(self.result_files):
            #if f.endswith('.resplus.json'):
                #datF.append(f.replace('.resplus.json',''))
        #datF.sort()
        #self.bp_lview_res.clear()
        #self.bp_lview_res.addItems(datF)

    def bp_buildMetaproteome(self):
        items = self.bp_lview.selectedItems()
        if len(items) == 0:
            self.bp_lview.selectAll()
            items = self.bp_lview.selectedItems()
        items2 = []
        for i in items:
            sp_str = str(i.text()).split(')-(')
            items2.append('{}'.format(sp_str[0][1:]))
        #print(items)
        #print(items2)

        #self.optimized_metaproteome = self.bp_btn_idopt.isChecked()
        self.optimized_metaproteome = self.enableOptimizationApp.isChecked()
        print('Optimization:', self.optimized_metaproteome)
        self.buildMetaProteomeFromSeqplus(items2, self.seqplus_files,
                                          self.optimized_metaproteome)
        r_html = self.buildHtmlStringMetaprot()
        self.widgetDisplayReact_update(r_html)
        self.bp_btn_blast.setEnabled(True)
        self.bp_btn_blast.update()

    def runBLASTFUNC(self):
        RUN_SEARCH = True
        input_fasta = str(self.bp_text_targ.text())
        if input_fasta == '' or input_fasta == None:
            print('\nERROR: Please load an input file!')
            return
        self.widget_busy.setLabelText("Running sequence search")
        self.widget_busy.setWindowTitle("Running sequence search")
        self.widget_busy.show()
        self.widgetBusyUpdate(10)
        wDir = os.path.join(self.blast_work_dir, str(time.time()).split('.')[0])
        linkd = self.link_file
        metap = self.metaproteome_file
        self.widgetBusyUpdate(20)
        inp_exec = os.path.join(self.blast_tools, 'inparanoid41_base.zip')
        #if os.sys.platform == 'win32':
            #inp_exec = os.path.join(self.blast_tools, 'inparanoid41_base_win_x64.zip')
        #else:
            #inp_exec = os.path.join(self.blast_tools, 'inparanoid41_base_linux.zip')
        self.widgetBusyUpdate(30)
        outgroup = None
        if self.bp_text_out.isEnabled():
            outgroup = str(self.bp_text_out.text())
            #print(outgroup)
            if not os.path.exists(outgroup):
                self.widgetBusyUpdate(40)
                self.status_bar.showMessage('Error in benchmark file: \"{}\"'.format(outgroup))
                #self.bp_text_out.setText('')
                outgroup = None
                self.menu_enableBenchmark()
        self.widgetBusyUpdate(50)
        _DAT_LINK_DICT_, resraw = self.runSequenceSearch(input_fasta, linkd, metap, wDir, inp_exec, outgroup, RUN_SEARCH)
        self.widgetBusyUpdate(80)
        if self.__DEL_BLAST_TMP__:
            shutil.rmtree(wDir, ignore_errors=True)
        self.widgetBusyUpdate(100)
        if _DAT_LINK_DICT_ is None and resraw is None:
            title = "About MetaDraft"
            msg = "There has been an error in the BLAST search subsystems\nplease see data_blast directory for debug information\n"
            self.widgetMsgBox(QMessageBox.Critical, title, msg)

    def widgetMsgBox(self, status, title, msg):
        """
        QMessageBox.Information
        QMessageBox.Critical
        QMessageBox.Warning
        """
        print('MSG: ' + msg)
        self._msg_box_ = QMessageBox(status, title, msg)
        self._msg_box_.show()


    def runInparanoid(self, par):
        wdir = par[0]
        target = par[1]
        dbase = par[2]
        out = par[3]
        end_flag = os.path.join(wdir, '.done')
        error_flag = os.path.join(wdir, '.fail')
        if os.path.exists(end_flag):
            os.remove(end_flag)

        os.chdir(wdir)
        if out is None:
            os_call = ['perl', 'inparanoid.pl', target, dbase]
        else:
            os_call = ['perl', 'inparanoid.pl', target, dbase, out]

        if self.__DEBUG__: print('\nWork directory: {}\nOS call: {}'.format(wdir, ' '.join(os_call)))

        #if raw_input('\nProceed with inParanoid? [y/n]: ')  == 'y':
        if True:
            TSTART = time.time()
            #self.blast_process = QtCore.QProcess()
            try:
                out = subprocess.check_call(os_call, shell=False)
            except subprocess.CalledProcessError as err:
                out = err.returncode
                if err.returncode == 2:
                    print('\n\nPERL/BLAST ERROR: possible no homology between source and target proteomes!')
                F = open(error_flag, 'w')
                F.close()
            #out = self.blast_process.start(os_call)
            TEND = time.time()

        runtime = int(math.floor((TEND-TSTART)/60.0))

        F = open(end_flag, 'w')
        F.close()

        print('\n\ninParanoid run took {} minutes to complete with return code: {}'.format(runtime, out))


    def runSequenceSearch(self, input_fasta, linkd, metap, wDir, inp_exec, outgroup, run_inparanoid):
        # if needed parse input genbank file and output a inParanoid fasta
        if input_fasta.endswith('.gbk') or input_fasta.endswith('.gb') or input_fasta.endswith('.gbff'):
            print('INFO: GenBank input detected creating FASTA file:')
            input_fasta = biotools.createParanoidFASTAfromFile([input_fasta])
            #print(input_fasta)

        # set up the inParanoid directory and input/output files
        psetup = setupInparanoid(inp_exec, wDir, input_fasta, metap, outgroup)
        jF = open(linkd, 'r')
        linkDict = json.load(jF)
        jF.close()

        linkDict['__metaproteome__']['input_fasta'] = input_fasta
        linkDict['__metaproteome__']['selection_state'] = {}
        #linkDict['__metaproteome__']['selection_state']['selected_genes'] = None
        #linkDict['__metaproteome__']['selection_state']['selected_reactions'] = None
        #linkDict['__metaproteome__']['selection_state']['selected_metabolites'] = None

        # run inParanoid
        #print(psetup)
        INP_START = time.time()
        if run_inparanoid:
            self.runInparanoid(psetup)
        INP_END = time.time()
        if os.path.exists(os.path.join(psetup[0], '.fail')):
            os.chdir(self.cDir)
            return None, None

        # parse inParanoid input/output
        inPtab = open(os.path.join(wDir, 'table.IN-DB'), 'r')

        input_seq_length = {}
        input_fasta_ids = []
        for seq_record in biotools.SeqIO.parse(input_fasta, "fasta"):
            #print(seq_record)
            input_seq_length[seq_record.id] = len(seq_record.seq)
            input_fasta_ids.append(seq_record.id)

        #setup results dictionary
        resraw = {}
        for l in inPtab:
            if l.startswith('OrtoID'):
                pass
            else:
                r = [a.strip() for a in l.split('\t')]
                r[2] = [a.strip() for a in r[2].split()][0]
                r[3] = [a.strip() for a in r[3].split()]
                GO = True
                pairs = []
                while GO:
                    a = r[3].pop(0)
                    b = float(r[3].pop(0))
                    pairs.append((a,b))
                    if len(r[3]) < 2:
                        r[3] = pairs
                        GO = False
                #print(r)

                resraw[r[2]] = {'id' : int(r[0]),
                                'bits' : int(r[1]),
                                'source' : r[2],
                                'match' : r[3],
                                #'length' : linkDict['__metaproteome__']['protein_lengths'][r[2]],
                                'length' : input_seq_length[r[2]],
                                'total' : linkDict['__metaproteome__']['total_length']}
        inPtab.close()

        resmatch = {}
        for r_ in resraw:
            resmatch[r_] = {}
            for m_ in resraw[r_]['match']:
                resmatch[r_][m_[0]] = m_[1]
        for gid in input_fasta_ids:
            if gid not in resmatch:
                resmatch[gid] = None
                print('INFO: match not found: {}'.format(gid))
                linkDict['__metaproteome__']['reports']['genes']['unmatched'].append(str(gid))

        linkDict['__metaproteome__']['search_results'] = resmatch
        if self.bp_btn_out.isEnabled():
            linkDict['__metaproteome__']['benchmark'] = str(self.bp_text_out.text())
        else:
            linkDict['__metaproteome__']['benchmark'] = None

        outfname1 = os.path.split(input_fasta)[-1]
        Fj = open(os.path.join(wDir, outfname1.replace('.fasta','.search_results.json')), 'w')
        json.dump(resmatch, Fj, indent=1, separators=(',', ': '))
        Fj.close()

        #Fj = open(linkd, 'w')
        new_ldict = os.path.join(wDir, os.path.split(linkd)[-1])
        Fj = open(new_ldict, 'w')
        json.dump(linkDict, Fj, indent=1, separators=(',', ': '))
        Fj.close()

        respath = self.result_files
        respath = os.path.join(respath, self.func_getCurrentUser(), time.strftime('%y-%m-%d'))
        if not os.path.exists(respath):
            os.makedirs(respath)
        outfname2 = os.path.split(linkd)[-1].replace('.json', '.resplus.json')
        if linkDict['__metaproteome__']['optimization']:
            outfname2 = outfname2.replace('_metalink', '-(opt)_metalink')
        if linkDict['__metaproteome__']['benchmark'] is not None:
            outfname2 = outfname2.replace('_metalink', '-(f)_metalink')
        outfname2 = outfname2.replace('.in.fasta','')
        #print(outfname1, outfname2)
        outfname1 = outfname1.replace('.in.fasta','')
        # update name to current result file name
        linkDict['__metaproteome__']['file_name'] = os.path.join(respath, '({})-{}'.format(outfname1, outfname2))
        Fj = open(os.path.join(respath, '({})-{}'.format(outfname1, outfname2)), 'w')
        json.dump(linkDict, Fj, indent=1, separators=(',', ': '))
        Fj.close()
        #writeLogBLAST(os.path.join(self.cDir, 'blast.log'), int(math.floor((INP_END - INP_START)/60.0)),\
                            #input_fasta, linkd, metap)
        print(input_fasta)
        os.remove(input_fasta)
        os.chdir(self.cDir)
        self.bp_text_targ.setText('')
        return linkDict, resraw

    def buildHtmlStringMetaprot(self):
        r_html = "<html><body>"
        r_html += "<table cellpadding=\"5\", border=\"1\", cellspacing=\"0\"><caption></caption>"
        r_html += "<tr><td>Input</td><td>{}</td></tr>".format(str(self.bp_text_targ.text()))
        r_html += "<tr><td>Benchmark</td><td>{}</td></tr>".format(str(self.bp_text_out.text()))
        r_html += "<tr><td>Metaproteome</td><td>{}</td></tr>".format(self.metaproteome_file)
        r_html += "<tr><td>Optimization</td><td>{}</td></tr>".format(self.optimized_metaproteome)
        r_html += "</table></body></html>"
        return r_html


    def buildMetaProteomeFromSeqplus(self, oid_list, model_lib, optimized_metaproteome):
        """
        Build a metaproteome from a list of organisms ids.

        - *oid_list* a list of one or more organism ids that connect to a seqplus model and link dictionary
        - *model_lib* the directory that contains the model lib

        """
        self.widget_busy.setWindowTitle("Building Metaproteome")
        self.widget_busy.setLabelText("Building Metaproteome")
        self.widget_busy.setMinimumWidth(100)
        self._progress_ = 10
        self.widgetBusyUpdate(self._progress_)
        self.widget_busy.show()

        outDir = self.metaproteome_files
        linkDict = self.buildLinkedDictFromSeqplus(oid_list, model_lib)
        #print(linkDict.keys())
        #print(oid_list)

        for o in oid_list:
            #old
            sbml1 = linkDict[o]["sbml_out"]
            sbml2 = linkDict[o]["sbml_out_generic"]

            if os.sys.platform == 'win32':
                p, f = os.path.split(sbml1)
                p, f = os.path.split(sbml2)
                linkDict[o]["sbml_out"] = os.path.join(self.seqplus_files, f)
                linkDict[o]["sbml_out_generic"] = os.path.join(self.seqplus_files, f)
            else:
                sbml1 = os.path.normpath(sbml1)
                sbml2 = os.path.normpath(sbml2)
                sbml1 = sbml1.split('\\')[-1]
                sbml2 = sbml2.split('\\')[-1]
                linkDict[o]["sbml_out"] = os.path.join(self.seqplus_files, sbml1)
                linkDict[o]["sbml_out_generic"] = os.path.join(self.seqplus_files, sbml2)

            #print(linkDict[o]["sbml_out"])
            #print(linkDict[o]["sbml_out_generic"])

        # "phylogenetic filtering" or "user defined ranking of organisims using reaction id's"
        if optimized_metaproteome:
            biotools.idFilter(linkDict, oid_list)

        oidString = '-'.join(oid_list)

        self._progress_ = 50
        self.widgetBusyUpdate(self._progress_)


        # create new style metaproteome
        fullname = oidString
        if len(oidString) > 100:
            oidString = 'metaproteome-{}'.format(time.strftime('%y%m%d%H%M'))
            print('Truncating metaproteome name to: {}').format(oidString)

        biotools.createMetaProteome(os.path.join(outDir, '({})_metaproteome.fasta'.format(oidString)),\
                                    linkDict, optimized=optimized_metaproteome, paranoid_style=True)
        linkDict['__metaproteome__']['optimization'] = optimized_metaproteome
        linkDict['__metaproteome__']['__fullname__'] = fullname

        self._progress_ = 90
        self.widgetBusyUpdate(self._progress_)

        if optimized_metaproteome:
            lfname = '({})_metalink.json'.format(oidString)
        else:
            lfname = '({})_metalink.json'.format(oidString)

        Fj = open(os.path.join(outDir, lfname), 'w')
        json.dump(linkDict, Fj, indent=1, separators=(',', ': '))
        Fj.close()

        #for n in range(len(biotools.IDMAP0)):
            #print(len(biotools.IDMAP0[n]), len(biotools.IDMAP1[n]))

        self.metaproteome = linkDict
        self.metaproteome_file = os.path.join(outDir, '({})_metaproteome.fasta'.format(oidString))
        self.link_file = os.path.join(outDir, '({})_metalink.json'.format(oidString))
        self._progress_ = 100
        self.widgetBusyUpdate(self._progress_)
        self.status_bar.showMessage('Metaproteome: \"{}\" created.'.format(self.metaproteome_file))

    def buildLinkedDictFromSeqplus(self, oid_list, model_lib):
        """
        Build a linked dict from a list of organisms ids and seqplus files.

        - *oid_list* a list of one or more organism ids that connect to a seqplus model and link dictionary
        - *model_lib* the directory that contains the model lib

        """
        assert os.path.exists(model_lib), '\nModel library not found at location: \"{}\"'.format(model_lib)
        files = os.listdir(model_lib)
        assert len(files) >= 2, '\nModel library must contain at least two files.'
        #print(files)
        linkDict = {}
        self._progress_ = 10
        for o_ in oid_list:
            for f_ in files:
                #print(o_)
                #print(f_.split(')-(', 1)[0][1:])
                if f_.split(')-(', 1)[0][1:] == o_ and f_.endswith('.json'):
                    #print(f_.split(')-(', 1)[0][1:])
                    F = open(os.path.join(model_lib, f_), 'r')
                    odict = json.load(F)
                    linkDict[o_] = odict[o_]
                    if '__idx__' in linkDict:
                        linkDict['__idx__'].update(odict['__idx__'])
                    else:
                        linkDict['__idx__'] = odict['__idx__']

                    F.close()
                    del odict
                    break
                self._progress_ += 5
                self.widgetBusyUpdate(self._progress_)

        return linkDict

    def widgetTableGene(self, parent=None):
        self.table_gene = QTableWidget(len(self._DAT_SEARCH_RES), len(self.table_gene_cols), parent=self.widget_tabpanel)
        self.table_gene.setColumnWidth(0, 100)
        self.table_gene.setHorizontalHeaderLabels(self.table_gene_cols)
        self.widgetTableGene_populate()
        self.table_gene.setSortingEnabled(True)
        self.table_gene.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.table_gene.customContextMenuRequested.connect(self.widgetTableGene_tableRightClicked)

        # this works for mouseclicks
        #self.table_gene.cellClicked.connect(self.widgetTableGene_cellClicked)
        # this should also work for cursor changes
        self.table_gene.itemSelectionChanged.connect(self.widgetTableGene_cellSelectionChanged)

        self.table_gene.cellChanged.connect(self.widgetTableGene_cellChecked)
        self._gene_all_ids_ = tuple([str(self.table_gene.item(r, 1).text()) for r in range(self.table_gene.rowCount())])
        self._gene_all_ids_src_ = tuple([str(self.table_gene.item(r, 0).text())[len(self.gene_prefix):] for r in range(self.table_gene.rowCount())])
        self._gene_all_ids_org_ = tuple([str(self.table_gene.item(r, 4).text()) for r in range(self.table_gene.rowCount())])
        self._gene_status_changed_ = False

    @pyqtSlot(QtCore.QPoint)
    def widgetTableGene_tableRightClicked(self, QPos):
        print('GENE TABLE RIGHT CLICKED', QPos)

        @pyqtSlot()
        def filterSelection():
            self.widget_tablegene_selectApp = QWidget()
            self.widget_tablegene_selectApp.setWindowTitle('Filter selected genes')
            ppos = self.func_getNewPopupWindowCoords()
            self.widget_tablegene_selectApp.setGeometry(QtCore.QRect(ppos[0], ppos[1], 350, 100))

            low_label = QLabel(self.widget_tablegene_selectApp)
            low_label.setText('Minimum score')
            up_label = QLabel(self.widget_tablegene_selectApp)
            up_label.setText('Maximum score')
            low_field = QLineEdit(self.widget_tablegene_selectApp)
            low_field.setText(str(self._gene_score_limits_[0]))
            up_field = QLineEdit(self.widget_tablegene_selectApp)
            up_field.setText(str(self._gene_score_limits_[1]))

            @pyqtSlot()
            def applySelection():
                input_error = False
                try:
                    low = float(str(low_field.text()))
                except ValueError:
                    low_field.setText('ERROR')
                    input_error = True
                try:
                    high = float(str(up_field.text()))
                except:
                    up_field.setText('ERROR')
                    input_error = True
                if not input_error:
                    if low > high:
                        input_error = True
                        self.widgetMsgBox(QMessageBox.Critical, 'Input error', 'low <= high')
                        return
                    elif low > 1.0 or high > 1.0:
                        input_error = True
                        self.widgetMsgBox(QMessageBox.Critical, 'Input error', 'low <= 1.0 and high <= 1.0')
                        return

                if input_error:
                    self.widgetMsgBox(QMessageBox.Critical, 'Input error', 'Invalid value')
                    return

                self._gene_score_limits_[0] = low
                self._gene_score_limits_[1] = high

                for gidx in range(self.table_gene.rowCount()):
                    try:
                        score = float(str(self.table_gene.item(gidx, 2).text()).strip())
                    except ValueError:
                        continue
                    checked = QtCore.Qt.Unchecked
                    if score >= low and score <= high:
                        checked = QtCore.Qt.Checked
                    self.table_gene.item(gidx, 3).setCheckState(checked)
                self.widget_tablegene_selectApp.close()


            applybut = QPushButton(self.widget_tablegene_selectApp)
            applybut.setText('Apply')
            applybut.clicked.connect(applySelection)

            layout = QGridLayout(self.widget_tablegene_selectApp)
            layout.setSpacing(10)
            layout.addWidget(low_label, 0, 0)
            layout.addWidget(up_label, 0, 1)
            layout.addWidget(low_field, 1, 0)
            layout.addWidget(up_field, 1, 1)
            layout.addWidget(applybut, 2, 0, 1, 2)

            self.widget_tablegene_selectApp.show()


        self.widget_tablegene_rclickmenu = QMenu()
        menu_item = self.widget_tablegene_rclickmenu.addAction("Filter selection")
        menu_item.triggered.connect(filterSelection)
        parentPosition = self.table_gene.viewport().mapToGlobal(QPos)
        self.widget_tablegene_rclickmenu.move(parentPosition)
        self.widget_tablegene_rclickmenu.show()

    def widgetTableGene_populate(self):
        row = 0
        sres = self._DAT_SEARCH_RES
        self._SRC_TRG_MAP_ = []
        # if needed in future
        ##self.table_gene.setRowCount(0)
        ##self.table_gene.setColumnCount(0)
        # maybe ...
        #self.table_gene.setRowCount(len(self._DAT_SEARCH_RES))
        #self.table_gene.setColumnCount(len(self.table_gene_cols))
        self.gene_table_loading = True
        for t_ in sres:
            if sres[t_] != None:
                grp_colour = None
                #print(t_, sres[t_])
                if len(sres[t_]) > 1:
                    grp_colour = next(self._colourCycler_)
                    vals = [sres[t_][g] for g in sres[t_]]
                else:
                    vals = [list(sres[t_].values())[0]]

                for g_ in sres[t_]:
                    item0 = QTableWidgetItem(str(t_))
                    item1 = QTableWidgetItem(str(g_))
                    item2 = NumberTableWidgetItem(str(sres[t_][g_]))
                    item3 = QTableWidgetItem()
                    item3.g_src = str(t_)
                    item3.g_match = str(g_)
                    item3.org = None
                    if g_ in self._DAT_LINK_DICT_['__idx__']:
                        org = str(self._DAT_LINK_DICT_['__idx__'][g_])
                        item4 = QTableWidgetItem(org)
                        item3.org = org
                    else:
                        item4 = QTableWidgetItem('NO GPR')
                    if item1 not in self._SRC_TRG_MAP_:
                        self._SRC_TRG_MAP_.append('{} {}'.format(item0, item1))

                    items = [item0, item1, item2, item3, item4]
                    for c_ in range(len(items)):
                        if c_ == 3:
                            items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                            items[c_].setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                            if sres[t_][g_] == max(vals):
                                items[c_].setCheckState(QtCore.Qt.Checked)
                            else:
                                items[c_].setCheckState(QtCore.Qt.Unchecked)
                        else:
                            items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                        if grp_colour != None:
                            items[c_].setBackground(QColor(*grp_colour))
                        self.table_gene.setItem(row, c_, items[c_])

                    row += 1
            else:
                self.table_gene.setItem(row, 0, QTableWidgetItem(str(t_)))
                self.table_gene.setItem(row, 1, QTableWidgetItem(''))
                self.table_gene.setItem(row, 2, NumberTableWidgetItem(''))
                chkBox = QTableWidgetItem()
                chkBox.setCheckState(QtCore.Qt.Unchecked)
                chkBox.setFlags(chkBox.flags() & ~QtCore.Qt.ItemIsUserCheckable & ~QtCore.Qt.ItemIsEnabled)
                self.table_gene.setItem(row, 3, chkBox)
                self.table_gene.setItem(row, 4, QTableWidgetItem(''))
                row += 1
        self.gene_table_loading = False
        self._updateGeneMap_()
        self._gene_status_changed_ = False
        self.table_gene.sortByColumn(0, Qt.SortOrder(0))

    def widgetTableMetab(self):
        self.table_metab = QTableWidget(0, len(self.table_metab_cols), parent=self.widget_tabpanel)
        self.table_metab.setColumnWidth(0,100)
        self.widgetTableMetab_populate()
        self.table_metab.setSortingEnabled(True)
        #self.table_metab.currentCellChanged.connect(self.widgetTableMetab_cellClicked)
        self.table_metab.cellClicked.connect(self.widgetTableMetab_cellClicked)
        self._metab_selected_map_ = self.widgetTableMetab_getMap()
        self._metab_status_changed_ = False

    def getMetabFromReactions(self):
        selected_metabolites = {}
        #print(self.widgetTableReaction_getMap())
        #print(self.widgetTableMetab_getMap())

        for r in self._reaction_selected_map_:
            if self._reaction_selected_map_[r]:
                R = self.selected_reactions[r]['obj']
                for S in R.getSpeciesObj():
                    #S = S.clone()
                    S._organism_ = R._organism_
                    if S not in selected_metabolites:
                        selected_metabolites[S.getId()] = S
        print('M', len(selected_metabolites))
        return selected_metabolites

    def widgetTableMetab_populate(self):
        self.table_metab.clear()
        self.table_metab.setRowCount(0)
        self.table_metab.setColumnCount(0)

        selected_metabolites = self.getMetabFromReactions()

        self.table_metab.setRowCount(len(selected_metabolites))
        self.table_metab.setColumnCount(len(self.table_metab_cols))
        self.table_metab.setHorizontalHeaderLabels(self.table_metab_cols)
        row = 0
        for s_ in selected_metabolites:
            S = selected_metabolites[s_]
            item0 = QTableWidgetItem(S.getId())
            item1 = QTableWidgetItem(S.getName())
            item2 = QTableWidgetItem(S._organism_)
            item3 = QTableWidgetItem()
            items = [item0, item1, item2, item3]
            for c_ in range(len(items)):
                if c_ == 3:
                    items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                    items[c_].setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                    if s_ not in self.selected_metabolites:
                        items[c_].setCheckState(QtCore.Qt.Unchecked)
                    elif s_ in self._metab_selected_map_:
                        if self._metab_selected_map_[s_]:
                            items[c_].setCheckState(QtCore.Qt.Checked)
                        else:
                            items[c_].setCheckState(QtCore.Qt.Unchecked)
                    else:
                        items[c_].setCheckState(QtCore.Qt.Unchecked)
                else:
                    items[c_].setFlags(items[c_].flags() & ~Qt.ItemIsEditable)
                self.table_metab.setItem(row, c_, items[c_])
            row += 1

        self.selected_metabolites = selected_metabolites
        #if len(self._metab_selected_map_) > 0:
        self._updateMetaboliteMap_()
        self._metab_status_changed_ = False
        self.table_metab.sortByColumn(0, Qt.SortOrder(0))

    @pyqtSlot(int, int)
    def widgetTableMetab_cellClicked(self, row, column):
        print("Metab: Row %d and Column %d was clicked" % (row, column))
        if self.table_metab.item(row, 0) is None:
            return
        self.METAB_LAST_SELECTED_ROW = row
        org = str(self.table_metab.item(row, 2).text())
        sid = str(self.table_metab.item(row, 0).text())
        metab = self.selected_metabolites[sid]
        r_html = self.buildHtmlStringsMetab(metab.getId())
        self.widgetDisplayReact_update(r_html)

        if self._wnotes_current_obj_ is not None:
            self.note_readFromNotesWidget(self._wnotes_current_obj_)
        self.note_writeToNotesWidget(metab, self.func_getCurrentUser(realname=False))
        self._wnotes_current_obj_ = metab

    def widgetTableMetab_getMap(self):
        mmap = {}
        for r in range(self.table_metab.rowCount()):
            mid = str(self.table_metab.item(r, 0).text())
            if self.table_metab.item(r, 3).checkState() == Qt.Checked:
                mmap[mid] = True
                self.selected_metabolites[mid].setBoundary()
            else:
                mmap[mid] = False
                self.selected_metabolites[mid].unsetBoundary()
        return mmap

    def widgetTableGene_getMap(self):
        mmap = {}
        for r in range(self.table_gene.rowCount()):
            mid = str(self.table_gene.item(r, 0).text())
            tid = str(self.table_gene.item(r, 1).text())
            key = '{}{}{}'.format(mid, self.id_sep, tid)
            if self.table_gene.item(r, 3).checkState() == Qt.Checked:
                mmap[key] = True
            else:
                mmap[key] = False
        return mmap

    def widgetTableGene_getSelectedIds(self):
        #tuple([str(self.table_gene.item(r, 1).text()) for r in range(self.table_gene.rowCount()) if self.table_gene.item(r, 3).checkState() == Qt.Checked])
        if self._gene_selected_map_ is not None:
            return [i.split(self.id_sep)[1] for i in self._gene_selected_map_ if self._gene_selected_map_[i]]
        else:
            return []

    def widgetTableReaction_getSelectedIds(self):
        #tuple([str(self.table_reaction.item(r, 0).text()) for r in range(self.table_reaction.rowCount()) if self.table_reaction.item(r, 2).checkState() == Qt.Checked])
        if self._reaction_selected_map_ is not None:
            return [i for i in self._reaction_selected_map_ if self._reaction_selected_map_[i]]
        else:
            return []

    def widgetTableReaction_getMap(self):
        rmap = {}
        for r in range(self.table_reaction.rowCount()):
            if self.table_reaction.item(r, 0) is not None and self.table_reaction.item(r, 2) is not None:
                rid = str(self.table_reaction.item(r, 0).text())
                org = str(self.table_reaction.item(r, 3).text())
                rkey = '{}{}{}'.format(rid, self.id_sep, org)
                if self.table_reaction.item(r, 2).checkState() == Qt.Checked:
                    rmap[rkey] = True
                else:
                    rmap[rkey] = False
        return rmap

    @pyqtSlot()
    def widgetTableGene_cellSelectionChanged(self):
        items = self.table_gene.selectedItems()
        #print(','.join([i.text() for i in items]))
        #print(items)
        if len(items) == 0:
            print('No items doing nothing')
        elif items[0].row() == self.GENE_LAST_SELECTED_ROW:
            print('Same row doing nothing')
        else:
            print('New row {} --> {}'.format(self.GENE_LAST_SELECTED_ROW, items[0].row()))
            self.widgetTableGene_cellClicked(items[0].row(), items[0].column())

    @pyqtSlot(int, int)
    def widgetTableGene_cellClicked(self, row, column):
        realmatch = True
        if self.table_gene.item(row, 0) == None:
            realmatch = False
        print("Gene: Row %d and Column %d was clicked" % (row, column))
        self.GENE_LAST_SELECTED_ROW = row
        gene = str(self.table_gene.item(row, 0).text())
        r_html = self.buildHtmlStringsGene(str(self.table_gene.item(row, 1).text()), gene)
        self.widgetDisplayReact_update(r_html)

        if realmatch:
            if self._wnotes_current_obj_ is not None:
                self.note_readFromNotesWidget(self._wnotes_current_obj_)
            gene = gene[len(self.gene_prefix):]
            self._dummy_gene_obj_.id = gene
            self._dummy_gene_obj_.setName(gene)
            self._dummy_gene_obj_.__setObjRef__(self._dummy_gene_obj_)
            self._wnotes_current_obj_ = self._dummy_gene_obj_
            self.note_writeToNotesWidget(self._dummy_gene_obj_, self.func_getCurrentUser(realname=False))

    @pyqtSlot(int, int)
    def widgetTableGene_cellChecked(self, row, column):
        if self.gene_table_loading:
            return
        #print("Gene: Row %d and Column %d was checked" % (row, column))
        itm = self.table_gene.item(row, column)
        src = str(self.table_gene.item(row, 0).text())
        gen = str(self.table_gene.item(row, 1).text())
        #print('GENEKEY: {}{}{}'.format(src, self.id_sep, gen))
        self.GENE_LAST_CHECKED_ROW = row
        self._updateGeneMap_(from_click=True)
        self.GENE_SELECTION_STATE_CHANGE = True
        #print(self._gene_selected_ids_)

    @pyqtSlot()
    def widgetTableReaction_cellSelectionChanged(self):
        items = self.table_reaction.selectedItems()
        #print(','.join([i.text() for i in items]))
        #print(items)
        if len(items) == 0:
            print('No items doing nothing')
        elif items[0].row() == self.REACT_LAST_SELECTED_ROW:
            print('Same row doing nothing')
        else:
            print('New row {} --> {}'.format(self.REACT_LAST_SELECTED_ROW, items[0].row()))
            self.widgetTableReaction_cellClicked(items[0].row(), items[0].column())

    @pyqtSlot(int, int)
    def widgetTableReaction_cellClicked(self, row, column):
        print("Reaction: Row %d and Column %d was clicked" % (row, column))
        if self.table_reaction.item(row, 0) == None:
            return
        self.REACT_LAST_SELECTED_ROW = row
        org = str(self.table_reaction.item(row, 3).text())
        reac = str(self.table_reaction.item(row, 0).text())
        reac = self._DAT_MODELS[org].getReaction(reac)
        r_html = self.buildHtmlStringsReaction(reac)
        self.widgetDisplayReact_update(r_html)
        if self._wnotes_current_obj_ is not None:
            self.note_readFromNotesWidget(self._wnotes_current_obj_)
        self.note_writeToNotesWidget(reac, self.func_getCurrentUser(realname=False))
        self._wnotes_current_obj_ = reac

    @pyqtSlot(int, int)
    def widgetTableReaction_cellChecked(self, row, column):
        if self.reaction_table_loading:
            return
        #print("Reaction: Row %d and Column %d was checked" % (row, column))
        itm = self.table_reaction.item(row, column)
        self.REACT_LAST_CHECKED_ROW = row
        self.REACTION_SELECTION_STATE_CHANGE = True
        self._updateReactionMap_(from_click=True)

    def getGeneAnnotationFromGeneDB(self, gid):
        annot = self._genedb_.getCell('GENES', 'id', gid, 'annotation')
        try:
            annot = json.loads(annot)
        except ValueError:
            annot = {}
        return annot

    def insertNotesToNotesDB(self, model, sid, notes, session):
        """
        self._notesdb_sqlcols_ = ['unixtime REAL PRIMARY KEY', 'user TEXT', 'model TEXT', 'id TEXT', 'notes TEXT', 'other TEXT']

        """
        data = {'unixtime' : time.time(),
                'user' : self.func_getCurrentUser(realname=False),
                'model' : str(model),
                'id' : str(sid),
                'notes' : str(notes),
                'other' : str(session)
                }
        #pprint.pprint(data)
        self._notesdb_.insertData('NOTES', data, commit=True)
        #if str(notes) != self._notesdb_.getCell('NOTES', 'notes', str(notes), 'notes'):
            #pprint.pprint(data)
            #self._notesdb_.insertData('NOTES', data, commit=True)

    def readNotesFromNotesDB(self, model, sid, user, session):
        sql = 'SELECT notes FROM NOTES WHERE model=\"{}\" AND id=\"{}\" AND user=\"{}\" AND other=\"{}\" ORDER BY unixtime DESC'.format(model, sid, user, session)
        notes = ''
        try:
            notes = self._notesdb_.fetchAll(sql)
            if len(notes) > 0:
                notes = [str(a[0]).strip() for a in notes][0]
            else:
                notes = ''
        except ValueError:
            notes = ''
        return notes

    def note_readFromNotesWidget(self, obj):
        #print('note_readFromNotesWidget')
        if obj is not None:
            notes = str(self._wnotes_.toPlainText())
            #print('note_readFromNotesWidget', notes)
            if notes not in ['', ' ']:
                #print('-inserting')
                #print(obj.__objref__().getName(), obj.getId())
                self.insertNotesToNotesDB(obj.__objref__().getName(), obj.getId(), notes, self._DAT_NOTESDB_KEY_)
            else:
                pass
                #print('-skipping empty1')
        else:
            pass
            #print('-skipping empty2')

    def note_writeToNotesWidget(self, obj, user):
        #print('note_writeToNotesWidget')
        if obj is not None:
            #print(obj.__objref__().getName(), obj.getId())
            notes = self.readNotesFromNotesDB(obj.__objref__().getName(), obj.getId(), user, self._DAT_NOTESDB_KEY_)
            #print('note_writeToNotesWidget', notes)
            self._wnotes_.setText(notes)

    def getCellX(self, table, col, rid, cell):
        """
        Get the table cell which correspond to rid in column. Returns the value or None

         - *table* the database table
         - *col* the column id
         - *rid* the row index id
         - *cell* the column of the cell you want tp extract

        """

        #print(sql)
        data = None
        try:
            data = str(self.db_cursor.execute(sql).fetchone()[0])
        except AttributeError:
            return None
        return data

    def getDBXrefFromGeneDB(self, gid):
        try:
            annot = self._genedb_.getCell('GENES', 'id', gid, 'db_xref')
            annot = self._genedb_.URLDecode(annot)
        except (ValueError, AttributeError):
            annot = {}
        return annot

    def func_formatGeneAnnotationToHTML(self, gene, gannot, r_html, color='#ffffff'):
        for ga_ in gannot:
            if ga_ in ['GO_process', 'GO_function', 'GO_component']:
                #match = re.search(self.regex['GOterm'], gannot[ga_])
                match = re.findall(self.regex['GOterm'], gannot[ga_])
                gstring = gannot[ga_]
                for m in range(len(match)):
                    href = '<a href=\"http://identifiers.org/go/{}\">{}</a>'.format(match[m], match[m])
                    if m >= 1:
                        href = '<br/>' + href
                    gstring = gstring.replace(match[m], href)
                r_html += "<tr bgcolor=\"{}\"><td><strong>{}</strong></td><td>{}</td></tr>".format(color, ga_, gstring)
            elif ga_ == 'EC_number':
                href = '<a href=\"http://identifiers.org/ec-code/{}\">{}</a>'.format(gannot[ga_], gannot[ga_])
                r_html += "<tr bgcolor=\"{}\"><td><strong>{}</strong></td><td>{}</td></tr>".format(color, ga_, href)
            elif ga_ in ['raw_location', 'codon_start', 'transl_table', 'locus_tag']:
                pass
            else:
                r_html += "<tr bgcolor=\"{}\"><td><strong>{}</strong></td><td>{}</td></tr>".format(color, ga_, gannot[ga_])
        try:
            db_xref =  eval(self.getDBXrefFromGeneDB(gene))
        except (TypeError, AttributeError):
            #print('dbxref error')
            db_xref = {}
        db_xrefs = ''
        for i in db_xref:
            if 'GeneID:' in i:
                db_xrefs += '<a href=\"http://www.ncbi.nlm.nih.gov/gene/{}\">{}</a><br/>'.format(i.split(':')[1], i)
            elif 'UniProtKB' in i:
                db_xrefs += '<a href=\"http://www.uniprot.org/uniprot/{}\">{}</a><br/>'.format(i.split(':')[1], i)
            elif 'GI:' in i:
                db_xrefs += '<a href=\"http://www.ncbi.nlm.nih.gov/protein/{}\">{}</a><br/>'.format(i.split(':')[1], i)
            else:
                db_xrefs += '{}<br/>'.format(i)
        r_html += "<tr bgcolor=\"{}\"><td><strong>{}</strong></td><td>{}</td></tr>".format(color, 'db_xref', db_xrefs)

        return r_html

    def buildHtmlStringsGene(self, gene, genesrc):
        print(gene, genesrc)
        realmatch = True
        if gene == '' or gene is None:
            realmatch = False
        genesrc = genesrc[len(self.gene_prefix):]
        if realmatch:
            organism = self._DAT_LINK_DICT_["__idx__"][gene]
            reactions = self._DAT_G2REACT[gene]
        #print(self._DAT_SEARCH_RES)
        #print(self._DAT_G2REACT)
        #self._updateGeneMap_()

        r_html = "<table cellpadding=\"5\", border=\"1\", cellspacing=\"0\"><caption>GeneMatch: {} --> {}</caption>".format(genesrc, gene)
        gannot = gannotsrc = None
        try:
            gannot = self.getGeneAnnotationFromGeneDB(gene)
        except TypeError:
            #print('gannot error')
            gannot = {}
        try:
            gannotsrc = self.getGeneAnnotationFromGeneDB(genesrc)
        except TypeError:
            #print('gannot error')
            gannotsrc = {}

        #print('\nDebug:', gene, genesrc)
        #print(gannot)
        #print(gannotsrc)

        if len(gannotsrc) > 0:
            r_html += '<tr><td colspan="2" align=\"center\"><strong>Source Gene Information: {}</strong></td></tr>'.format(genesrc)
            r_html = self.func_formatGeneAnnotationToHTML(genesrc, gannotsrc, r_html, color='#e6ffe6')


        if realmatch:
            r_html += '<tr><td colspan="2" align=\"center\"><strong>Matching Gene Information: {}</strong></td></tr>'.format(gene)
            r_html = self.func_formatGeneAnnotationToHTML(gene, gannot, r_html, color='#ffffcc')
            r_html += '<tr><td colspan="2" align=\"center\"><strong>Associated reactions(s)</strong></td></tr>'

            for o_ in reactions:
                #r_html += "<tr><td><strong>Organism</strong></td><td><strong>{}</strong></td></tr>".format(organism)
                #if o_.reversible:
                    #reverse = '(reversible)'
                #else:
                    #reverse = '(irreversble)'
                #r_html += "<tr><td><strong>{}</strong><br/>{}</td><td><strong>{}</strong></td></tr>".format(o_.getId(), reverse, o_.getName())
                r_html += "<tr><td><strong>{}</strong> ({})</td><td><strong>{}</strong></td></tr>".format(o_.getId(), organism, o_.getName())
                #r_html += "<tr><td><strong>Equation ID</strong></td><td>{}</td></tr>".format(o_.getEquation())
                #r_html += "<tr><td><strong>Equation Name</strong></td><td>{}</td></tr>".format(o_.getEquation(use_names=True))
                #r_html += "<tr><td><strong>Source gene</strong></td><td>{}</td></tr>".format(gene)
                #r_html += "<tr><td><strong>Matching gene</strong></td><td>{}</td></tr>".format(gene)
                try:
                    emod = self._DAT_MODELS[o_._organism_]
                    gpr = emod.getGPRforReaction(o_.getId())
                    if emod.__FBC_VERSION__ < 2:
                        assoc = gpr.getAssociationStr(use_labels=True)
                    else:
                        assoc = gpr.getAssociationStr(use_labels=True)
                    # TODO: use a regular expression here
                    gids = gpr.getGeneLabels()
                    gids = list(set(gids))
                    gids.sort(key = len)
                    gids.reverse()
                    #print(self._gene_selected_ids_)
                    #print(self._gene_selected_map_)
                    #print(assoc)
                    #print(gids)
                    if len(gids) == 1:
                        if gids[0] in self._gene_selected_ids_:
                            assoc = assoc.replace(assoc, '<span style="color: green;">{}</span> '.format(gids[0]))
                        elif gids[0] in self._gene_all_ids_:
                            assoc = assoc.replace(assoc, '<span style="color: red;">{}</span> '.format(gids[0]))
                    else:
                        for g_ in gids:
                            #ids = str('{}{}{}'.format(genesrc, self.id_sep, g_))
                            if g_ in self._gene_selected_ids_:
                                #print(ids, g_ in self._gene_selected_ids_)
                                assoc = assoc.replace(g_+' ', '<span style="color: green;">{}</span> '.format(g_))
                                assoc = assoc.replace(' '+g_, ' <span style="color: green;">{}</span>'.format(g_))
                                #assoc = assoc.replace(g_, ' <span style="color: green;">{}</span>'.format(g_))
                            elif g_ in self._gene_all_ids_:
                                assoc = assoc.replace(g_+' ', '<span style="color: red;">{}</span> '.format(g_))
                                assoc = assoc.replace(' '+g_, ' <span style="color: red;">{}</span>'.format(g_))
                                #assoc = assoc.replace(g_, ' <span style="color: red;">{}</span>'.format(g_))
                    assoc = assoc.replace(gene, '<strong>{}</strong>'.format(gene))

                except AttributeError as why:
                    o_.serializeToDisk(o_.getId())
                    print('ERROR: {} - {}'.format(gene, o_.getId()))
                    print(why)
                    assoc = '<span style="color: red;"><strong>{}</strong></span> '.format('UNKNOWN')

                #r_html += "<tr><td><strong>Association</strong></td><td>{}</td></tr>".format(assoc)
                r_html += '<tr><td colspan="2" align=\"left\">{}</td></tr>'.format(assoc)
                # disabled for now, using geneDB
                #miriam = o_.getMIRIAMannotations()
                #if miriam != None:
                    #for m in miriam:
                        #if len(miriam[m]) > 0:
                            #for u in range(len(miriam[m])):
                                #r_html += "<tr><td>{}</td><td><a href=\"{}\">{}</a></td></tr>".format(m, miriam[m][u], miriam[m][u])
        r_html += "</table>"

        return r_html

    def buildHtmlStringsReaction_getGeneInfo(self, gene):
        html = ''
        try:
            gannot = self.getGeneAnnotationFromGeneDB(gene)
        except TypeError:
            gannot = {}
        for ga_ in gannot:
            html += "<strong>{}: </strong>{}<br/>".format(ga_, gannot[ga_])
        return html

    def buildHtmlStringsReaction(self, reac):
        if reac == '':
            return 'No matching reaction'

        bkg_clr_sub = '#CCFFCC'
        bkg_clr_prd = '#FFCCCC'

        r_html = '<html><head>'

        r_html += '</head><body>'
        r_html += "<table><tr><td><strong>Organism</strong></td><td><strong>{}</strong></td></tr>".format(reac._organism_)
        if reac.reversible:
            reverse = '(reversible)'
        else:
            reverse = '(irreversble)'
        r_html += "<tr><td><strong>{}</strong><br/>{}</td><td><strong>{}</strong></td></tr>".format(reac.getId(), reverse, reac.getName())
        _gene_selected_ids_ = self.widgetTableGene_getSelectedIds()
        gpr = self._DAT_MODELS[reac._organism_].getGPRforReaction(reac.getId())
        assoc = gpr.getAssociationStr(use_labels=True)
        assoc_new = ''
        #tuple([str(self.table_reaction.item(r, 0).text()) for r in range(self.table_reaction.rowCount()) if self.table_reaction.item(r, 2).checkState() == Qt.Checked])
        for r in range(self.table_reaction.rowCount()):
            rkey = '{}{}{}'.format(str(self.table_reaction.item(r, 0).text()), self.id_sep, str(self.table_reaction.item(r, 3).text()))
            if rkey == '{}{}{}'.format(reac.getId(), self.id_sep, reac._organism_):
                assoc_new = str(self.table_reaction.item(r, 5).text())
                if ',' in assoc_new:
                    # this might be or
                    assoc_new = ' and '.join(list(set(assoc_new.split(','))))
                assoc_new = '({})'.format(assoc_new)
                #print(assoc_new)
                break

        # TODO: use a regular expression here
        gids = gpr.getGeneLabels()
        gids.sort()
        gids.reverse()

        for g_ in gids:
            ginfo = self.buildHtmlStringsReaction_getGeneInfo(g_)
            if g_ in _gene_selected_ids_:
                ginfo = self.buildHtmlStringsReaction_getGeneInfo(g_)
                if len(gids) == 1:
                    assoc = assoc.replace(g_, '<span title="{}" style="color: green;">{}</span> '.format(ginfo, g_))
                else:
                    assoc = assoc.replace(g_+' ', '<span title="{}" style="color: green;">{}</span> '.format(ginfo, g_))
                    assoc = assoc.replace(' '+g_, ' <span title="{}" style="color: green;">{}</span>'.format(ginfo, g_))
                    #assoc = '<p title="{}">{}</p>'.format(g_, assoc)
            elif g_ in self._gene_all_ids_:
                ginfo = self.buildHtmlStringsReaction_getGeneInfo(g_)
                assoc = assoc.replace(g_+' ', '<span title="{}" style="color: red;">{}</span> '.format(ginfo, g_))
                assoc = assoc.replace(' '+g_, ' <span title="{}" style="color: red;">{}</span>'.format(ginfo, g_))
        r_html += "<tr><td><strong>Association</strong></td><td>{}</td></tr>".format(assoc)
        r_html += "<tr><td><strong>Association new</strong></td><td>{}</td></tr>".format(assoc_new)
        r_html += "<tr><td><strong>Equation ID</strong></td><td>{}</td></tr>".format(reac.getEquation())
        r_html += "<tr><td><strong>Equation Name</strong></td><td>{}</td></tr>".format(reac.getEquation(use_names=True))
        r_html += '<tr><td colspan="2" align=\"center\"><strong>Substrates</strong></td></tr>'
        for s in reac.getSubstrateIds():
            S = reac.__objref__().getSpecies(s)
            r_html += "<tr><td bgcolor=\"{}\"><strong>{}</strong></td><td bgcolor=\"{}\"><strong>{}</strong></td></tr>".format(bkg_clr_sub, S.getId(), bkg_clr_sub, S.getName())
            r_html += "<tr><td>charge: {}</td><td>{}</td></tr>".format(S.getCharge(), S.getChemFormula())
        r_html += '<tr><td colspan="2" align=\"center\"><strong>Products</strong></td></tr>'
        for p in reac.getProductIds():
            P = reac.__objref__().getSpecies(p)
            r_html += "<tr><td bgcolor=\"{}\"><strong>{}</strong></td><td bgcolor=\"{}\"><strong>{}</strong></td></tr>".format(bkg_clr_prd, P.getId(), bkg_clr_prd, P.getName())
            r_html += "<tr><td>charge: {}</td><td>{}</td></tr>".format(P.getCharge(), P.getChemFormula())
        r_annot = reac.getAnnotations()
        r_html += '<tr><td colspan="2" align=\"center\"><strong>Annotations</strong></td></tr>'
        for ant in r_annot:
            if not ant.startswith('gbank_seq_'):
                r_html += "<tr><td>{}:</td><td>{}</td></tr>".format(ant, r_annot[ant])
        r_html += "</table>"
        r_html += '</body></html>'
        return r_html

    def buildHtmlStringsMetab(self, metab):
        if metab == '':
            return 'No matching metabolite'
        bkg_clr_sub = '#FFFFFF'
        S = self.selected_metabolites[metab]
        r_html = '<html><body>'
        r_html += '<table cellpadding="5", border="1", cellspacing="0"><caption></caption>'
        r_html += "<tr><td bgcolor=\"{}\"><strong>{}</strong></td><td bgcolor=\"{}\"><strong>{}</strong></td></tr>".format(bkg_clr_sub, 'Id', bkg_clr_sub, S.getId())
        r_html += "<tr><td bgcolor=\"{}\"><strong>{}</strong></td><td bgcolor=\"{}\"><strong>{}</strong></td></tr>".format(bkg_clr_sub, 'Name', bkg_clr_sub, S.getName())
        r_html += "<tr><td>Charge</td><td>{}</td></tr>".format(S.getCharge())
        r_html += "<tr><td>Formula</td><td>{}</td></tr>".format(S.getChemFormula())
        reag_of = [r for r in S.isReagentOf() if r in [rr.split('@')[0] for rr in self.selected_reactions]]
        r_html += "<tr><td>ReagentOf</td><td>{}</td></tr>".format(', '.join(reag_of))
        #S.addMIRIAMannotation('is', 'CheBI', 'CHEBI:17822')
        #S.addMIRIAMannotation('isPartOf', 'CheBI', 'CHEBI:17822')
        #S.addMIRIAMannotation('is', 'CheBI', 'CHEBI:17822')
        miriam = S.getMIRIAMannotations()
        if miriam != None:
            r_html += "<tr><td colspan=\"2\" align=\"center\">RDF references (opens in browser)</td></tr>"
            for m in miriam:
                if len(miriam[m]) > 0:
                    for u in range(len(miriam[m])):
                        r_html += "<tr><td>{}</td><td><a href=\"{}\">{}</a></td></tr>".format(m, miriam[m][u], miriam[m][u])
        r_html += '</table></body></html>'
        return r_html

    def widgetResultTree(self):
        self.tree_results = QTreeWidget()
        self.tree_results.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tree_results.setUniformRowHeights(True)
        self.tree_results.setColumnCount(1)
        #self.tree_results.setColumnWidth(0, 300)
        self.tree_results.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tree_results.customContextMenuRequested.connect(self.widgetTreeRightClickMenu)
        self.tree_results.setHeaderLabels(['Results ({})'.format(self.func_getCurrentUser(realname=True))])
        self.tree_results.itemSelectionChanged.connect(self.bp_rtreeOnSelect)
        self.widgetResultTree_fill()

    @pyqtSlot(QtCore.QPoint)
    def widgetTreeRightClickMenu(self, position):
        indexes = self.tree_results.selectedIndexes()
        treeidx = self.tree_results.indexAt(position)
        if len(indexes) == 1:
            print(indexes[0].parent().data(Qt.DisplayRole))
            if indexes[0].parent().data(Qt.DisplayRole) is None:
                self._widget_result_tree_rightclick_data_ = str(indexes[0].data(Qt.DisplayRole))
            else:
                self._widget_result_tree_rightclick_data_ = os.path.join(str(indexes[0].parent().data(Qt.DisplayRole)),
                                                                         str(indexes[0].data(Qt.DisplayRole)))
        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1

        menu = QMenu(self.tree_results)
        if level == 0:
            menu_item0 = menu.addAction(self.tr("Rename"))
            menu_item0.triggered.connect(self.widgetResultTreeRightClickRename)
            menu_item1 = menu.addAction(self.tr("Delete"))
            menu_item1.triggered.connect(self.widgetResultTreeRightClickDelete)
        elif level == 1:
            menu_item1 = menu.addAction(self.tr("Delete"))
            menu_item1.triggered.connect(self.widgetResultTreeRightClickDelete)
        menu.exec_(self.tree_results.viewport().mapToGlobal(position))

    @pyqtSlot()
    def widgetResultTreeRightClickRename(self):
        path = os.path.join(self.result_files, self.func_getCurrentUser(), self._widget_result_tree_rightclick_data_)
        text, ok = QInputDialog.getText(self, 'Input Dialog',
                                              'Enter new directory name:')
        if ok:
            path2 = path.replace(self._widget_result_tree_rightclick_data_, text)
            print(path)
            print(path2)
            if not os.path.exists(path2):
                shutil.move(path,  path2)
                self.widgetResultTree_fill()
            else:
                self.widgetMsgBox(QMessageBox.Warning, 'Warning Dialog', 'Directory {} already exists'.format(text))

    @pyqtSlot()
    def widgetResultTreeRightClickDelete(self):
        path = os.path.join(self.result_files, self.func_getCurrentUser(), self._widget_result_tree_rightclick_data_)
        #print(os.path.split(self._widget_result_tree_rightclick_data_))
        DELTREE = True

        if os.path.split(self._widget_result_tree_rightclick_data_)[0] != '':
            path += '_metalink.resplus.json'
            DELTREE = False

        if os.path.exists(path):
            reply = QMessageBox.question(self, 'Message',\
                                               "Are you sure you want to delete:\n{}?".format(self._widget_result_tree_rightclick_data_),\
                                               QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                if not DELTREE:
                    os.remove(path)
                else:
                    shutil.rmtree(path)
                self.widgetResultTree_fill()
                print('Deleted file \"{}\"'.format(path))

    def widgetResultTree_item(self, item, value):
        if type(value) is dict:
            for key, val in value.items():
                child = QTreeWidgetItem()
                if val is None:
                    child.setText(0, u'{}'.format(key.replace('_metalink.resplus.json', '')))
                    ## this is the timestamp code that has now been removed
                    ##print(item.text(0), key, val)
                    ##try:
                        ##path = os.path.join(self.getRtreeItemPath(item), key)
                        ##t = datetime.datetime.fromtimestamp(os.path.getmtime(path))
                        ##child.setText(1, unicode('{}:{}'.format(t.hour, t.minute)))
                    ##except WindowsError:
                        ##print('DEBUG: {}'.format(self.getRtreeItemPath(item)), key, value)
                    child.setExpanded(True)
                else:
                    child.setText(0, u'{}'.format(key))

                item.addChild(child)
                self.widgetResultTree_item(child, val)

    def bp_rtreeOnSelect(self):
        getSelected = self.tree_results.selectedItems()
        if getSelected:
            path = self.getRtreeItemPath(getSelected[0])
            path += '_metalink.resplus.json'
            #print(path)
            if os.path.exists(path):
                self.result_file = path
                #print('Good path {}'.format(path))
            else:
                self.result_file = ''
                #print('Invalid path {}'.format(path))


    def getRtreeItemPath(self, item):
        path = os.path.join(self.result_files, self.func_getCurrentUser())
        item_path = self.func_getTreePathOfItem(item)
        for i_ in item_path:
            path = os.path.join(path, i_)
        return path


    def func_getTreePathOfItem(self, item):
        path = []
        while item is not None:
            try:
                path.insert(0, str(item.text(0)))
                item = item.parent()
            except:
                print('func_getTreePathOfItem error')
                item = None
        return path

    def func_getCurrentUser(self, realname=False):
        if not realname:
            return self._CONFIG_['system']['current_user']
        else:
            return self._CONFIG_['users'][self._CONFIG_['system']['current_user']]['name']

    def widgetResultTree_fill(self):
        self.tree_results.clear()
        #print(self.result_files)
        root = os.path.join(self.result_files, self.func_getCurrentUser())
        root = self.getDirectoryStructure(root)[self.func_getCurrentUser()]
        #print(root)
        self.widgetResultTree_item(self.tree_results.invisibleRootItem(), root)
        self.tree_results.expandToDepth(2)

    def getDirectoryStructure(self, rootdir):
        """
        Creates a nested dictionary that represents the folder structure of rootdir
        """
        dirout = {}
        rootdir = rootdir.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in os.walk(rootdir):
            folders = path[start:].split(os.sep)
            subdir = dict.fromkeys(files)
            parent = reduce(dict.get, folders[:-1], dirout)
            parent[folders[-1]] = subdir
        return dirout

    def widgetTabPanel(self):
        self.widget_tabpanel = QTabWidget(parent=self)
        self.widget_tabpanel.currentChanged.connect(self.onTabChange)

    def widgetTabPanel_createtab(self, widget):
        tab = QWidget()
        tab.layout = QVBoxLayout()
        tab.layout.addWidget(widget)
        tab.setLayout(tab.layout)
        return tab

    def widgetTabPanel_add(self, widget, tabname, setcurrent=False):
        self.widget_tabpanel.addTab(widget, tabname)
        if setcurrent:
            self.widget_tabpanel.setCurrentWidget(widget)

    def widgetTabPanelRight(self):
        self.widget_tabpanel_right = QTabWidget(self)
        self.widget_tabpanel_right.currentChanged.connect(self.onTabChange)

    def widgetTabPanelRight_createtab(self, widget):
        tab = QWidget()
        tab.layout = QVBoxLayout()
        tab.layout.addWidget(widget)
        tab.setLayout(tab.layout)
        return tab

    def widgetTabPanelRight_add(self, widget, tabname, setcurrent=False):
        self.widget_tabpanel_right.addTab(widget, tabname)
        if setcurrent:
            self.widget_tabpanel_right.setCurrentWidget(widget)

    def _update_Reactions_(self):
        self.widgetTableReaction_populate()
        self._gene_status_changed_ = False
        #gselect = self.widgetTableGene_getSelectedIds()
        #if self._gene_selected_ids_ != gselect:
            #self._gene_status_changed_ = True
            #self._gene_selected_ids_ = gselect
            #self.widgetTableReaction_populate()
        #else:
            #self._gene_status_changed_ = False

    def _update_Metabolites_(self):
        self.widgetTableMetab_populate()


    def _updateReactionMap_(self, from_click=False):
        self._reaction_selected_map_ = self.widgetTableReaction_getMap()
        self._reaction_selected_ids_ = self.widgetTableReaction_getSelectedIds()
        if from_click:
            row = self.REACT_LAST_CHECKED_ROW
            self.status_bar.showMessage('Selected Reactions: {}'.format(len(self._reaction_selected_ids_)))

            #gene = str(self.table_reaction.item(row, 1).text())
            #genesrc = str(self.table_reaction.item(row, 0).text())
            #print(row, gene, genesrc)

        #print('R', len(self._reaction_selected_map_), len(self._reaction_selected_ids_))
        #print(self._reaction_selected_map_)
        #print(self._reaction_selected_ids_)

    def _updateMetaboliteMap_(self):
        self._metab_selected_map_ = self.widgetTableMetab_getMap()

    def _updateGeneMap_(self, from_click=False):
        self._gene_selected_map_ = self.widgetTableGene_getMap()
        self._gene_selected_ids_ = self.widgetTableGene_getSelectedIds()
        if from_click:
            row = self.GENE_LAST_CHECKED_ROW
            gene = str(self.table_gene.item(row, 1).text())
            genesrc = str(self.table_gene.item(row, 0).text())
            #print(row, gene, genesrc)
            html = self.buildHtmlStringsGene(gene, genesrc)
            self.widgetDisplayReact_update(html)
            self.status_bar.showMessage('Selected Genes: {}'.format(len(self._gene_selected_ids_)))
            #print(self._gene_selected_map_)
            #print(self._gene_selected_ids_)
            #print('G', len(self._gene_selected_map_), len(self._gene_selected_ids_))

    @pyqtSlot(int)
    def onTabChange(self, idx):
        if self._loading_:
            return

        QApplication.setOverrideCursor(QCursor(QtCore.Qt.WaitCursor))
        self._last_tab_ = self._active_tab_
        self._active_tab_ = idx
        #print('tabChange', self._last_tab_, self._active_tab_)

        # read notes to object and clear text
        self.note_readFromNotesWidget(self._wnotes_current_obj_)
        self._wnotes_current_obj_ = None
        self._wnotes_.setText('')

        old_sel_gene = self._gene_selected_ids_
        old_sel_reac = self._reaction_selected_ids_

        if self._tabpanel_idx_[self._last_tab_] == 'Reactions':
            self._updateReactionMap_()
        elif self._tabpanel_idx_[self._last_tab_] == 'Metabolites':
            self._updateMetaboliteMap_()
        elif self._tabpanel_idx_[self._last_tab_] == 'Genes':
            print('Leaving genes')
            self._updateGeneMap_()
        elif self._tabpanel_idx_[self._last_tab_] == 'Build':
            pass

        update_generic_status = False

        if self.GENE_SELECTION_STATE_CHANGE:
            print('Updating reactions + metabolites')
            self._update_Reactions_()
            self._update_Metabolites_()
            self.GENE_SELECTION_STATE_CHANGE = False
            self.REACTION_SELECTION_STATE_CHANGE = False
        elif self.REACTION_SELECTION_STATE_CHANGE:
            print('Updating metabolites')
            self._update_Metabolites_()
            self.REACTION_SELECTION_STATE_CHANGE = False
        else:
            print('OnTabChange no selection state change')


        if self._tabpanel_idx_[self._active_tab_] == 'Genes':
            update_generic_status = True
            self._wnotes_.setEnabled(True)
            self.table_gene.clearSelection()
            self.table_gene.clearFocus()
            print('gene', self.GENE_LAST_SELECTED_ROW)
            self.table_gene.setCurrentCell(self.GENE_LAST_SELECTED_ROW, 0)
        elif self._tabpanel_idx_[self._active_tab_] == 'Reactions':
            print('Landing on reactions')
            update_generic_status = True
            self._wnotes_.setEnabled(True)
            #self.table_reaction.clearSelection()
            #self.table_reaction.clearFocus()
            print('react', self.REACT_LAST_SELECTED_ROW)
            self.table_reaction.setCurrentCell(self.REACT_LAST_SELECTED_ROW, 0)
        elif self._tabpanel_idx_[self._active_tab_] == 'Metabolites':
            self.table_metab.update()
            update_generic_status = True
            self._wnotes_.setEnabled(True)
            #self.status_bar.showMessage('Including {} metabolites.'.format(self.table_metab.rowCount()))
            print('metab', self.METAB_LAST_SELECTED_ROW)
            self.table_metab.setCurrentCell(self.METAB_LAST_SELECTED_ROW, 0)

        if update_generic_status:
            self.status_bar.showMessage('Genes: {} Reactions: {} Metabolites: {}'.format(len(self._gene_selected_ids_),\
                                                                                         len(self._reaction_selected_ids_), self.table_metab.rowCount()))

        self.func_saveSelectionState(False, False, False)
        QApplication.restoreOverrideCursor()

    def func_saveSelectionState(self, write, new_session, update_tables=True):
        if '__tmp__' not in self._DAT_LINK_DICT_['__metaproteome__']['selection_state']:
            self._DAT_LINK_DICT_['__metaproteome__']['selection_state']['__tmp__'] = {}

        if update_tables:
            self._updateGeneMap_()
            self._update_Reactions_()
            self._update_Metabolites_()

        selstate = self._DAT_LINK_DICT_['__metaproteome__']['selection_state']['__tmp__'].copy()
        selstate['selected_genes'] = self._gene_selected_map_
        selstate['selected_reactions'] = self._reaction_selected_map_
        selstate['selected_metabolites'] = self._metab_selected_map_
        if new_session:
            key, ok = QInputDialog.getText(self, 'Input Dialog',
                                                  'Enter session name:')
            key = str(key)
            if not ok or key in ['', ' ']:
                key = time.strftime('%Y-%m-%d-%H-%M-%S')
            self._DAT_LINK_DICT_['__metaproteome__']['selection_state'][key] = selstate
        if write:
            self.func_saveResultsFile()
        #print('I AM SAVING STATE')
        #pprint.pprint(selstate['selected_genes'])

    def func_saveResultsFile(self):
        fp = open(self.result_file, 'w')
        json.dump(self._DAT_LINK_DICT_, fp, indent=1, separators=(',', ': '))
        fp.flush()
        fp.close()

    def func_loadSelectionState(self, state):
        state = str(state)
        #print(state, self._DAT_LINK_DICT_['__metaproteome__']['selection_state'].keys())
        if state in self._DAT_LINK_DICT_['__metaproteome__']['selection_state']:
            selstate = self._DAT_LINK_DICT_['__metaproteome__']['selection_state'][state]
            self._gene_selected_map_ = copy.deepcopy(selstate['selected_genes'])
            self._reaction_selected_map_ = copy.deepcopy(selstate['selected_reactions'])
            self._metab_selected_map_ = copy.deepcopy(selstate['selected_metabolites'])
            print('I HAVE SET THE SELECTION STATE: {}'.format(state))
            self.CURRENT_SELECTION_STATE = state
            return True
        else:
            return False

    def func_setSelectionState(self):
        # set gene selected state
        print('self.table_gene.rowCount()', self.table_gene.rowCount())
        for row in range(self.table_gene.rowCount()):
            mid = str(self.table_gene.item(row, 0).text())
            tid = str(self.table_gene.item(row, 1).text())
            #print(mid,tid)
            key = '{}{}{}'.format(mid, self.id_sep, tid)
            #print(key, self._gene_selected_map_[key])
            if key not in self._DAT_LINK_DICT_['__metaproteome__']['selection_state'][self.CURRENT_SELECTION_STATE]['selected_genes']:
                print('StateSetError: {}'.format(key))
            else:
                if self._DAT_LINK_DICT_['__metaproteome__']['selection_state'][self.CURRENT_SELECTION_STATE]['selected_genes'][key]:
                    self.table_gene.item(row, 3).setCheckState(QtCore.Qt.Checked)
                else:
                    self.table_gene.item(row, 3).setCheckState(QtCore.Qt.Unchecked)
        self.table_gene.update()
        self.menu_buildAll()

        #self.widgetTableReaction_populate()
        ## set reaction selected state
        #self.widgetTableMetab_populate()
        ## set metab selection state
        #print('I HAVE SET THE SELECTION STATE')


    @pyqtSlot(int)
    def onTabRightChange(self, idx):
        self._last_tab_right_ = self._active_tab_right_
        self._active_tab_right_ = idx
        #print('r', self._last_tab_right_, self._active_tab_right_)
        #print(self._tabpanel_right_idx_[self._last_tab_right_], self._tabpanel_right_idx_[idx])

    def widgetPanel(self):
        self.widgetTableGene()
        self.widgetTableReaction()
        self.widgetTableMetab()

        tab1L = self.widgetTabPanel_createtab(self.table_gene)
        self.widgetTabPanel_add(tab1L, 'Genes', setcurrent=True)
        tab2L = self.widgetTabPanel_createtab(self.table_reaction)
        self.widgetTabPanel_add(tab2L, 'Reactions', setcurrent=False)
        tab3L = self.widgetTabPanel_createtab(self.table_metab)
        self.widgetTabPanel_add(tab3L, 'Metabolites', setcurrent=False)
        #tab4L = self.widgetTabPanel_createtab(QWidget())
        #self.widgetTabPanel_add(tab4L, 'Objective', setcurrent=False)
        #tab5L = self.widgetTabPanel_createtab(QWidget())
        #self.widgetTabPanel_add(tab5L, 'Other', setcurrent=False)

        # Dev
        #self.widget_tabpanel.setTabEnabled(0, True)
        self.widget_tabpanel.setCurrentIndex(1)

        self.widget_tabpanel.setTabEnabled(1, False)
        self.widget_tabpanel.setTabEnabled(2, False)
        self.widget_tabpanel.setTabEnabled(3, False)
        #self.widget_tabpanel.setTabEnabled(4, False)
        #self.widget_tabpanel.setTabEnabled(5, False)

        #self.grid.update()


class InputValidators(object):
    """Various input validators"""

    def inpv_stringItemNotInList(self, itm, items):
        if itm not in items:
            return True
        else:
            return False



    def inpv_floatItemInRange(self, itm, rmin, rmax):
        try:
            itm = float(itm)
            rmin = float(rmin)
            rmax = float(rmax)
        except:
            return False
        if rmin <= itm <= rmax:
            return True
        else:
            return False



class ConfigPanelWidgetINP(QWidget, InputValidators):
    """Generates a configurations panel from a dictionary of options"""

    kobjdict = None
    kdesc = None

    def __init__(self, dictobject, dictname, kdescript=None):
        super(ConfigPanelWidgetINP, self).__init__()
        self.kdesc = kdescript
        self.kobjdict = {}
        self.setWindowModality(Qt.ApplicationModal)

        keys =  list(getattr(dictobject, dictname).keys())
        keys.sort()

        self.grid = QGridLayout(self)
        self.grid.setSpacing(10)

        # config tooltips
        tooltips = {'PY_score_cutoff': 'Value must in range 0 <= value <= 100',
                    'PY_grey_zone': 'Value must in range 0 <= value <= 1',
                    'PY_conf_cutoff': 'Value must in range 0 <= value <= 1',
                    'PY_seq_overlap_cutoff': 'Value must in range 0 <= value <= 1',
                    'PY_group_overlap_cutoff': 'Value must in range 0 <= value <= 1',
                    'PY_segment_coverage_cutoff': 'Value must in range 0 <= value <= 1',
                    'PY_outgroup_cutoff': 'Value must in range 0 <= value <= 100',
                    'PY_matrix': 'One of: BLOSUM45, BLOSUM62, BLOSUM62, BLOSUM80, PAM70, PAM30'
                    }


        for r in range(len(keys)):
            if keys[r] not in ['PY_use_bootstrap', 'PY_use_outgroup']:
                k = QLabel(parent=self)
                k.setText(keys[r])
                k.setToolTip(tooltips[keys[r]])
                self.grid.addWidget(k, r, 0, 1, 1)

                v = QLineEdit(parent=self)
                v.setMaximumHeight(25)
                v.setText(getattr(dictobject, dictname)[keys[r]])
                v.mtk_keyid = keys[r]
                v.setToolTip(tooltips[keys[r]])
                self.kobjdict[keys[r]] = v
                self.grid.addWidget(v, r, 1, 1, 1)

        def bp_SaveExitFunc():
            pal = QPalette()
            textbad = QColor(Qt.red)
            textgood = QColor(Qt.black)

            GO = True
            for o in self.kobjdict:
                val = str(self.kobjdict[o].text())
                print(self.kobjdict[o].mtk_keyid, val)
                pal.setColor(QPalette.Text, textgood)
                self.kobjdict[o].setPalette(pal)
                if self.kobjdict[o].mtk_keyid == 'PY_matrix':
                    if not val in ["BLOSUM45","BLOSUM62","BLOSUM62","BLOSUM80","PAM70","PAM30"]:
                        GO = False
                        print('bad1', val)
                        pal.setColor(QPalette.Text, textbad)
                        self.kobjdict[o].setPalette(pal)
                # TODO: add type checks for numerical values
                elif self.kobjdict[o].mtk_keyid in ['PY_outgroup_cutoff', 'PY_score_cutoff']:
                    if not self.inpv_floatItemInRange(val, 0, 100):
                        GO = False
                        print('bad2', val)
                        pal.setColor(QPalette.Text, textbad)
                        self.kobjdict[o].setPalette(pal)
                elif self.kobjdict[o].mtk_keyid in ['PY_conf_cutoff', 'PY_seq_overlap_cutoff', 'PY_group_overlap_cutoff', 'PY_segment_coverage_cutoff', 'PY_grey_zone']:
                    if not self.inpv_floatItemInRange(val, 0, 1):
                        GO = False
                        print('bad3', val)
                        pal.setColor(QPalette.Text, textbad)
                        self.kobjdict[o].setPalette(pal)
                if GO:
                    getattr(dictobject, dictname)[o] = val
            self.update()
            if GO:
                self.close()

        bp_saveExit = QPushButton(parent=self)
        bp_saveExit.setText('Save and Exit')
        bp_saveExit.clicked.connect(bp_SaveExitFunc)
        self.grid.addWidget(bp_saveExit, len(keys)+1, 0, 1, 1)

        bp_Exit = QPushButton(parent=self)
        bp_Exit.setText('Exit')
        bp_Exit.clicked.connect(self.close)
        self.grid.addWidget(bp_Exit, len(keys)+1, 1, 1, 1)

        self.show()

BQBIOL = {'encodes' : 'encodes',
          'has part' : 'hasPart',
          'has version' : 'hasVersion',
          'has taxon' : 'hasTaxon',
          'is' : 'is',
          'is encoded by' : 'isEncodedBy',
          'is homolog to' : 'isHomologTo',
          'is part of' : 'isPartOf',
          'is version of' : 'isVersionOf'
          }

BQMODEL = {'is derived from' : 'isDerivedFrom',
           'is described by' : 'isDescribedBy',
           'occurs in' :  'occursIn'
           }

class CBMPyAnnotationEditor(QWidget):
    mod = None
    sid = None
    miriam = None
    bqbiol = None
    #bqmodel = None
    #bqmodelrev = None
    bqbiolrev = None
    fixed_font = QFont('Courier', 9.5)
    _fixColour = QColor(0,0,153,alpha=255)
    _errColour = QColor(255,0,0,alpha=255)
    _goodColour = QColor(0,100,0,alpha=255)
    _wnotes_ = None

    def __init__(self, mod, sid):
        super(CBMPyAnnotationEditor, self).__init__()
        self.setWindowModality(Qt.ApplicationModal)
        self.mod = mod
        self.sid = sid
        if sid in self.mod.getReactionIds():
            self.obj = self.mod.getReaction(sid)
        elif sid in self.mod.getSpeciesIds():
            self.obj = self.mod.getSpecies(sid)
        elif sid in self.mod.getGeneIds():
            self.obj = self.mod.getGene(sid)
        elif sid in self.mod.getCompartmentIds():
            self.obj = self.mod.getCompartment(sid)
        else:
            self.obj = cbmpy.CBModel.Fbase()
            self.obj.setId('id{}'.format(int(time.time())))

        self.miriam = copy.deepcopy(cbmpy.miriamids.miriamids)
        self.bqbiol = BQBIOL
        self.bqbiol.update(BQMODEL)
        #self.bqmodel = BQMODEL

        self.miriam_keys = list(self.miriam.keys())
        self.miriam_keys.sort()
        self.bqbiol_keys = list(self.bqbiol.keys())
        self.bqbiol_keys.sort()
        #self.bqmodel_keys = list(self.bqmodel.keys())
        #self.bqmodel_keys.sort()

        self.bqbiolrev = {value:key for key,value in self.bqbiol.items()}
        #self.bqmodelrev = {value:key for key,value in BQMODEL.items()}

        self.url2miriam = {}
        for k in self.miriam:
            self.url2miriam[self.miriam[k]['url']] = k

        self.initUI()

    @pyqtSlot()
    def btn_clicked_saveandclose(self):
        go = self.func_saveMiriam()
        if not go:
            return
        self.func_saveKVPairs()
        self.obj.setNotes(self.widget_lblnotes.toPlainText())
        self.close()

    def func_saveKVPairs(self):
        self.obj.annotation = {}
        for r in range(self.tblKvp.rowCount()):
            if not (str(self.tblKvp.item(r, 0).text()) == '' and str(self.tblKvp.item(r, 1).text()) == ''):
                self.obj.setAnnotation(str(self.tblKvp.item(r, 0).text()), str(self.tblKvp.item(r, 1).text()))

    def func_saveMiriam(self):
        newanno = []
        GO = True
        badanno = []
        for r in range(self.tblDesc.rowCount()):
            #try:
            qual = rsrc = idx = None
            if self.tblDesc.item(r, 0) is not None:
                qual = self.bqbiol[str(self.tblDesc.item(r, 0).text())]
                rsrc = str(self.tblDesc.item(r, 1).text())
                idx = str(self.tblDesc.item(r, 2).text())
            else:
                qual = self.bqbiol[str(self.tblDesc.cellWidget(r, 0).currentText())]
                rsrc = str(self.tblDesc.cellWidget(r, 1).currentText())
                idx = str(self.tblDesc.item(r, 2).text())
            if qual is not None and rsrc is not None and idx is not None:
                try:
                    if not self.obj.miriam.checkId(rsrc, idx):
                        GO = False
                        badanno.append((qual, rsrc, idx))
                    else:
                        newanno.append((qual, rsrc, idx))
                except:
                    GO = False
                    badanno.append((qual, rsrc, idx))


        #except Exception as e:
                #print(e)
                #print('Error with MIRIAM annotation row: {}'.format(r))
        if not GO:
            print(badanno)
            title = "MIRIAM identifier check"
            msg = 'ID check failures:\n'
            for fail in badanno:
                msg += '- \"{}\" is of class \"{}\" should have form: {}\n'.format(fail[2], fail[1], self.miriam[fail[1]]['example'])
            self._widget_msgBox_ = QMessageBox(QMessageBox.Warning, title, msg)
            self._widget_msgBox_.show()
            return False
        self.obj.miriam = None
        for qual, rsrc, idx in newanno:
            self.obj.addMIRIAMannotation(qual, rsrc, idx)
        return True


    def func_insertDescRow(self, qual, rsrc, sid, new=True):
        print('insertRow')
        newrow = self.tblDesc.rowCount()
        self.tblDesc.insertRow(newrow)
        if new:
            bqbiol = self.createBQBiolComboBoxItem()
            miriam = self.createMirimComboBoxItem()
            bqbiol.setCurrentIndex(self.bqbiol_keys.index(qual))
            miriam.setCurrentIndex(self.miriam_keys.index(rsrc))
            self.tblDesc.setCellWidget(newrow, 0, bqbiol)
            self.tblDesc.setCellWidget(newrow, 1, miriam)
        else:
            item0 = QTableWidgetItem(qual)
            item0.setFlags(QtCore.Qt.ItemIsSelectable |  QtCore.Qt.ItemIsEnabled )
            item1 = QTableWidgetItem(rsrc)
            item1.setFlags(QtCore.Qt.ItemIsSelectable |  QtCore.Qt.ItemIsEnabled )
            self.tblDesc.setItem(newrow, 0, item0)
            self.tblDesc.setItem(newrow, 1, item1)

        item2 = QTableWidgetItem(sid)
        self.tblDesc.setItem(newrow, 2, item2)

    def func_insertKvpRow(self, k, v):
        print('insertKVP')
        newrow = self.tblKvp.rowCount()
        self.tblKvp.insertRow(newrow)
        item0 = QTableWidgetItem(k)
        item1 = QTableWidgetItem(v)
        self.tblKvp.setItem(newrow, 0, item0)
        self.tblKvp.setItem(newrow, 1, item1)

    @pyqtSlot()
    def btn_clicked_add_desc(self):
        print('I CLICKED ADD DESC BUTTON')
        newrow = self.tblDesc.rowCount()
        self.tblDesc.insertRow(newrow)
        self.tblDesc.setCellWidget(newrow, 0, self.createBQBiolComboBoxItem())
        self.tblDesc.setCellWidget(newrow, 1, self.createMirimComboBoxItem())
        item = QTableWidgetItem('xxxxxxxxx')
        self.tblDesc.setItem(newrow, 2, item)
        self.tblDesc.resizeColumnsToContents()
        self.tblDesc.update()
        item.setText('')

    #@pyqtSlot()
    #def btn_clicked_add_create(self):
        #print('I CLICKED ADD CREATE BUTTON')

    #@pyqtSlot()
    #def btn_clicked_add_ref(self):
        #print('I CLICKED ADD REF BUTTON')

    @pyqtSlot()
    def btn_clicked_add_kvp(self):
        print('I CLICKED ADD KVP BUTTON')
        self.func_insertKvpRow('', '')


    def createMirimComboBoxItem(self):
        widget = QComboBox()
        widget.setMaximumWidth(200)
        cntr = 0
        for m_ in self.miriam_keys:
            widget.addItem('{}'.format(m_))
        return widget

    def createBQBiolComboBoxItem(self):
        widget = QComboBox()
        for m_ in self.bqbiol_keys:
            widget.addItem('{}'.format(m_))
        widget.view
        return widget

    #def createBQModelComboBoxItem(self):
        #widget = QComboBox()
        #for m_ in self.bqmodel_keys:
            #widget.addItem('{}'.format(m_))
        #return widget

    @pyqtSlot(QtCore.QPoint)
    def tblDesc_rightClickMenu(self, QPos):
        print('Description TABLE RIGHT CLICKED', QPos)
        # TODO try get a default table menu thing going

        try:
            idx = self.tblDesc.selectedIndexes()[0]
        except:
            return
        #self.tblDesc.selectRow(idx.row())
        val = str(self.tblDesc.item(idx.row(), 2).text())
        print(idx.row(), val)


        def deleteAnnotation():
            print('Delete annotation')
            self.tblDesc.removeRow(idx.row())

        def checkId():
            print('Check annotation identifier')
            entity = self.tblDesc.item(idx.row(), 1)
            if entity is None:
                entity = str(self.tblDesc.cellWidget(idx.row(), 1).currentText())
            else:
                entity = str(entity.text())
            title = "MIRIAM ID check"
            if self.obj.miriam is None:
                self.obj.miriam = cbmpy.CBDataStruct.MIRIAMannotation()

            if not self.obj.miriam.checkId(str(entity), val):
                print('Failed id check')
                msg = 'Invalid ID! Id should be of the form: {}'.format(self.miriam[entity]['example'])
                self._widget_msgBox_ = QMessageBox(QMessageBox.Warning, title, msg)
            else:
                msg = 'ID check passed'
                self._widget_msgBox_ = QMessageBox(QMessageBox.Information, title, msg)
            self._widget_msgBox_.show()

        self.widget_tbldesc_rclickmenu = QMenu()

        delr = self.widget_tbldesc_rclickmenu.addAction('Delete Annotation')
        delr.triggered.connect(deleteAnnotation)

        checkid = self.widget_tbldesc_rclickmenu.addAction('Check identifier')
        checkid.triggered.connect(checkId)

        #self.widget_tbldesc_rclickmenu.addSeparator()

        parentPosition = self.tblDesc.viewport().mapToGlobal(QPos)
        self.widget_tbldesc_rclickmenu.move(parentPosition)
        self.widget_tbldesc_rclickmenu.show()

    @pyqtSlot(QtCore.QPoint)
    def tblKvp_rightClickMenu(self, QPos):
        print('KVP TABLE RIGHT CLICKED', QPos)
        # TODO try get a default table menu thing going

        try:
            idx = self.tblKvp.selectedIndexes()[0]
        except:
            return
        #self.tblKvp.selectRow(idx.row())
        val = str(self.tblKvp.item(idx.row(), 1).text())
        print(idx.row(), val)


        def deleteAnnotation():
            print('Delete annotation')
            self.tblKvp.removeRow(idx.row())

        self.widget_tblkvp_rclickmenu = QMenu()

        delr = self.widget_tblkvp_rclickmenu.addAction('Delete Annotation')
        delr.triggered.connect(deleteAnnotation)

        parentPosition = self.tblKvp.viewport().mapToGlobal(QPos)
        self.widget_tblkvp_rclickmenu.move(parentPosition)
        self.widget_tblkvp_rclickmenu.show()

    def initUI(self):
        # create labels
        max_label_height = 30

        lbl_notes = QLabel(parent=self)
        lbl_notes.setText('Notes')

        # define label font from default font
        lbl_bld_font = lbl_notes.font()
        #lbl_bld_font.setBold(True)
        lbl_bld_font.setPointSize(10)

        lbl_notes.setFont(lbl_bld_font)
        lbl_notes.setAlignment(Qt.AlignLeft)
        lbl_notes.setMaximumHeight(max_label_height)

        lbl_desc = QLabel(parent=self)
        lbl_desc.setFont(lbl_bld_font)
        lbl_desc.setText('Description')
        lbl_desc.setAlignment(Qt.AlignLeft)
        lbl_desc.setMaximumHeight(max_label_height)

        #lbl_auth = QLabel(parent=self)
        #lbl_auth.setText('Creator')
        #lbl_auth.setFont(lbl_bld_font)
        #lbl_auth.setAlignment(Qt.AlignLeft)
        #lbl_auth.setMaximumHeight(max_label_height)

        #lbl_ref = QLabel(parent=self)
        #lbl_ref.setFont(lbl_bld_font)
        #lbl_ref.setText('References')
        #lbl_ref.setAlignment(Qt.AlignLeft)
        #lbl_ref.setMaximumHeight(max_label_height)

        lbl_kvp = QLabel(parent=self)
        lbl_kvp.setFont(lbl_bld_font)
        lbl_kvp.setText('CBMPy annotation')
        lbl_kvp.setAlignment(Qt.AlignLeft)
        lbl_kvp.setMaximumHeight(max_label_height)

        #lbl_sel = QLabel(parent=self)
        #lbl_sel.setFont(lbl_bld_font)
        #lbl_sel.setText('')
        #lbl_sel.setAlignment(Qt.AlignCenter)
        #lbl_sel.setMaximumHeight(max_label_height)

        # create tables
        table_min_width = 600
        self.tblDesc = QTableWidget(parent=self)
        try:
            self.tblDesc.horizontalHeader().ResizeMode(QtWidgets.QHeaderView.Stretch)
        except:
            self.tblDesc.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.tblDesc.horizontalHeader().setStretchLastSection(True)
        self.tblDesc.setMinimumWidth(table_min_width)
        #self.tblDesc.setFont(self.fixed_font)
        self.tblDesc.setShowGrid(True)
        self.tblDesc.setSortingEnabled(False)
        self.tblDesc.insertColumn(0)
        self.tblDesc.insertColumn(1)
        self.tblDesc.insertColumn(2)
        #self.tblDesc.insertColumn(3)
        self.tblDesc.setHorizontalHeaderLabels(['Relationship', 'Resource', 'ID'])
        self.tblDesc.verticalHeader().setVisible(True)
        self.tblDesc.horizontalHeader().setVisible(True)
        #self.tblDesc.itemClicked.connect(self.action_reactionProdSelect)
        #self.tblDesc.cellChanged.connect(self.checkBalance)
        self.tblDesc.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tblDesc.customContextMenuRequested.connect(self.tblDesc_rightClickMenu)

        #self.tblAuth = QTableWidget(parent=self)
        #self.tblAuth.setMinimumWidth(table_min_width)
        ##self.tblAuth.setFont(self.fixed_font)
        #self.tblAuth.setShowGrid(True)
        #self.tblAuth.setSortingEnabled(False)
        #self.tblAuth.insertColumn(0)
        #self.tblAuth.insertColumn(1)
        #self.tblAuth.insertColumn(2)
        #self.tblAuth.insertColumn(3)
        #self.tblAuth.setHorizontalHeaderLabels(['Famliy name', 'Given name', 'Email', 'Organisation'])
        #self.tblAuth.verticalHeader().setVisible(True)
        #self.tblAuth.horizontalHeader().setVisible(True)
        ##self.tblAuth.itemClicked.connect(self.action_reactionConSelect)
        #self.tblAuth.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        ##self.tblSub.cellChanged.connect(self.checkBalance)
        ##self.tblSub.customContextMenuRequested.connect(self.tblSub_rightClickMenu)
        #self.tblAuth.update()

        #self.tblRef = QTableWidget(parent=self)
        #self.tblRef.setMinimumWidth(table_min_width)
        ##self.tblRef.setFont(self.fixed_font)
        #self.tblRef.setShowGrid(True)
        #self.tblRef.setSortingEnabled(False)
        #self.tblRef.insertColumn(0)
        #self.tblRef.insertColumn(1)
        #self.tblRef.insertColumn(2)
        #self.tblRef.setHorizontalHeaderLabels(['Resource', 'ID', 'Description'])
        #self.tblRef.verticalHeader().setVisible(True)
        #self.tblRef.horizontalHeader().setVisible(True)
        ##self.tblRef.itemClicked.connect(self.action_reactionSpeciesSelect)
        #self.tblRef.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        #self.tblRef.update()

        table_min_width = 400
        self.tblKvp = QTableWidget(parent=self)
        self.tblKvp.setMinimumWidth(table_min_width)
        try:
            self.tblKvp.horizontalHeader().ResizeMode(QtWidgets.QHeaderView.Stretch)
        except:
            self.tblKvp.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        self.tblKvp.horizontalHeader().setStretchLastSection(True)
        #self.tblKvp.setFont(self.fixed_font)
        self.tblKvp.setShowGrid(True)
        self.tblKvp.setSortingEnabled(False)
        self.tblKvp.insertColumn(0)
        self.tblKvp.insertColumn(1)
        self.tblKvp.setHorizontalHeaderLabels(['Key', 'Value'])
        self.tblKvp.verticalHeader().setVisible(True)
        self.tblKvp.horizontalHeader().setVisible(True)
        #self.tblKvp.itemClicked.connect(self.action_reactionSpeciesSelect)
        self.tblKvp.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tblKvp.customContextMenuRequested.connect(self.tblKvp_rightClickMenu)
        self.tblKvp.update()

        ## populate tables
        if self.obj.miriam is not None:
            uris = self.obj.miriam.getAllMIRIAMUris()
            for u in uris:
                if len(uris[u]) > 0:
                    qual = self.bqbiolrev[u]
                    for uri in uris[u]:
                        rsrc = uri.split('/')
                        sid = rsrc[-1]
                        rsrc = self.url2miriam['/'.join(rsrc[:-1]) + '/']
                        self.func_insertDescRow(qual, rsrc, sid, new=False)
        self.tblDesc.resizeColumnsToContents()
        self.tblDesc.update()

        kvannot = self.obj.getAnnotations()
        if kvannot is not None:
            for k in kvannot:
                self.func_insertKvpRow(k, kvannot[k])
        self.tblKvp.resizeColumnToContents(0)
        self.tblKvp.update()

        ## Create buttons
        btn_saveandclose = QPushButton(parent=self)
        btn_saveandclose.setText('Save and exit')
        btn_saveandclose.clicked.connect(self.btn_clicked_saveandclose)

        btn_max_h = 25
        btn_max_w = 50

        btn_add_desc = QPushButton(parent=self)
        btn_add_desc.setText('Add')
        btn_add_desc.clicked.connect(self.btn_clicked_add_desc)
        btn_add_desc.setMaximumWidth(btn_max_w)
        btn_add_desc.setMaximumHeight(btn_max_h)

        #btn_add_ref = QPushButton(parent=self)
        #btn_add_ref.setText('Add')
        #btn_add_ref.clicked.connect(self.btn_clicked_add_ref)
        #btn_add_ref.setMaximumWidth(btn_max_w)
        #btn_add_ref.setMaximumHeight(btn_max_h)

        #btn_add_create = QPushButton(parent=self)
        #btn_add_create.setText('Add')
        #btn_add_create.clicked.connect(self.btn_clicked_add_create)
        #btn_add_create.setMaximumWidth(btn_max_w)
        #btn_add_create.setMaximumHeight(btn_max_h)

        btn_add_kvp = QPushButton(parent=self)
        btn_add_kvp.setText('Add')
        btn_add_kvp.clicked.connect(self.btn_clicked_add_kvp)
        btn_add_kvp.setMaximumWidth(btn_max_w)
        btn_add_kvp.setMaximumHeight(btn_max_h)

        ## create textFields
        self.widget_lblnotes = QTextEdit(parent=self)
        #self.widget_lblnotes.setMaximumHeight(80)
        #self.widget_lblnotes.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.widget_lblnotes.customContextMenuRequested.connect(self.modelDetail_RightClicked)
        self.widget_lblnotes.setHtml('{}'.format(self.obj.getNotes()))
        #self.widget_lblnotes.setReadOnly(True)
        self.widget_lblnotes.setAlignment(QtCore.Qt.AlignLeft)

        ## resize columns
        #self.tblReactProd.resizeColumnsToContents()
        #self.tblReactProd.resizeRowsToContents()
        #self.tblReactCon.resizeColumnsToContents()
        #self.tblReactCon.resizeRowsToContents()
        #self.tblSpecies.resizeColumnsToContents()
        #self.tblSpecies.resizeRowsToContents()

        # do layout
        #grid = QVBoxLayout(self)
        grid = QGridLayout(self)
        grid.setSpacing(10)

        #grid.addWidget(lbl_sel, 0, 0, 1, 2)
        grid.addWidget(lbl_desc, 0, 0, 1, 2)
        grid.addWidget(btn_add_desc, 0, 2, 1, 1)
        grid.addWidget(self.tblDesc, 1, 0, 11, 3)
        #grid.addWidget(lbl_ref, 4, 0, 1, 2)
        #grid.addWidget(btn_add_ref, 4, 2, 1, 1)
        #grid.addWidget(self.tblRef, 5, 0, 3, 3)
        #grid.addWidget(lbl_auth, 8, 0, 1, 2)
        #grid.addWidget(btn_add_create, 8, 2, 1, 1)
        #grid.addWidget(self.tblAuth, 9, 0, 3, 3)
        grid.addWidget(lbl_kvp, 0, 3, 1, 2)
        grid.addWidget(btn_add_kvp, 0, 5, 1, 1)
        grid.addWidget(self.tblKvp, 1, 3, 5, 3)
        grid.addWidget(lbl_notes, 6, 3, 1, 2)
        grid.addWidget(self.widget_lblnotes, 7, 3, 5, 3)
        grid.addWidget(btn_saveandclose, 12, 0, 1, 6)

        self.setLayout(grid)
        self.setWindowTitle('Editing annotation: {} - ({})'.format(self.sid, self.obj.getName()))
        self.show()



def setupInparanoid(ip_src, work_dir, target_fasta, metaproteome, outgroup):
    """
    This sets up an inParanoid session,
    """
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    assert os.path.exists(work_dir), "\nWork directory [{}] does not exist".format(work_dir)
    assert os.path.exists(ip_src), "\nCannot find inparanoid"

    assert os.path.exists(metaproteome), "Cannot find metaproteome [{}]".format(metaproteome)
    assert os.path.exists(target_fasta), "Cannot find target fasta [{}]".format(target_fasta)
    #assert os.path.exists(link), "Cannot find link dictionary [{}]".format(link)

    zfile = zipfile.ZipFile(ip_src, 'r')
    zfile.extractall(work_dir)
    zfile.close()
    inoid = os.path.join(work_dir, 'inparanoid.pl')
    if outgroup is None:
        bionoid.CONFIGKEYS['PY_use_outgroup'] = '0'
    else:
        bionoid.CONFIGKEYS['PY_use_outgroup'] = '1'
    bionoid.USERWIN, bionoid.USERLINUX = bionoid.buildUser(bionoid.CONFIGKEYS, bionoid.WINKEYS, bionoid.LINUXKEYS)

    if os.sys.platform == 'win32':
        bionoid.writeInparanoidBase(inoid, bionoid.HEAD, bionoid.USERWIN, bionoid.BODYWIN)
    else:
        bionoid.writeInparanoidBase(inoid, bionoid.HEAD, bionoid.USERLINUX, bionoid.BODYLINUX)

    para_in = 'IN'
    para_db = 'DB'
    para_out = None

    shutil.copyfile(target_fasta, os.path.join(work_dir, para_in))
    shutil.copyfile(metaproteome, os.path.join(work_dir, para_db))
    if outgroup is not None:
        para_out = 'TEST'
        shutil.copyfile(outgroup, os.path.join(work_dir, para_out))

    try:
        st = os.stat(os.path.join(work_dir, 'inparanoid.pl'))
        os.chmod(os.path.join(work_dir, 'inparanoid.pl'), st.st_mode | stat.S_IEXEC)
        os.chmod(os.path.join(work_dir, 'blast_parser.pl'), st.st_mode | stat.S_IEXEC)
    except:
        print('Could not change mode')

    return (work_dir, para_in, para_db, para_out)


def writeLogBLAST(log, inp_time, input_fasta, metalink, metap):
    if not os.path.exists(log):
        logF = open(log, 'w')
    else:
        logF = open(log, 'a')
    logF.write('blast,run,{},\"{}\",\"{}\",\"{}\"\n'.format(inp_time, input_fasta, metalink, metap))
    logF.close()


class MetaDraftApp(QMainWindow):
    def __init__(self):
        super(MetaDraftApp, self).__init__()
        self.initUI()

    def initUI(self):
        self._gui_ = MetaDraftGUI(self)
        self.setCentralWidget(self._gui_)
        self.setGeometry(200,200,1200,700)
        self.show()

def main():
    app = QApplication(sys.argv)
    widget_splash = QSplashScreen(QPixmap("images/metatoolkit1-03.jpg"))
    widget_splash.show()
    if RELEASE_STATUS == 2:
        widget_splash.showMessage("Ver {}-({}) beta\nAuthor: Brett G. Olivier PhD\n(c) Brett G. Olivier, Amsterdam, 2017.\nSee Help - About for more details\nb.g.olivier@vu.nl".format(__version__, cbmpy.__version__), alignment=Qt.AlignBottom)
        time.sleep(7)
    else:
        widget_splash.showMessage("Ver {}-({}) beta\n(c) Brett G. Olivier, Amsterdam, 2017.\nSee Help - About for more details\nb.g.olivier@vu.nl".format(__version__, cbmpy.__version__), alignment=Qt.AlignBottom)
        if RELEASE_STATUS == 1:
            time.sleep(2)
    ex = MetaDraftApp()
    widget_splash.finish(ex)
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

    #import profile
    #profile.run('main()', 'profile.tmp')

    #import pstats
    #p = pstats.Stats('profile.tmp')
    #p.sort_stats('cumulative').print_stats(50)
    #p.sort_stats('time').print_stats(50)

