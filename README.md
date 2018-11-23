# cbmpy-metadraft
CBMPy Metadraft: a flexible and extensible, GUI-based genome-scale model reconstruction tool that supports multiple Systems Biology standards.

## Installation and system dependencies

### System dependencies

MetaDraft requires working versions of Python, Perl, Java and NCBI Blast2. Once you have downloaded MetaDraft and it's Python dependencies, please run the *system test* script to that can be found in the MetaDraft directory `python systemtest.py`.

Typically the default Java runtime (JRE) installed on your system is all that is require, for Windows this can be obtained from [https://www.java.com/]. On Linux please use your favourite package manager, for example, when using Ubuntu/Debian based oeprating systems, try `sudo apt-get install default-jre`. 

MetaDraft has been successfully tested using Windows 10, Ubuntu Linux 16.04/18.04, Python 2.7, Python 3.6, PyQt4 and PyQt5. MetaDraft is developed on Windows Anaconda Python 2.7 with PyQt5.

### Getting MetaDraft

To install the latest version of MetaDraft, download from GitHub ([https://github.com/SystemsBioinformatics/cbmpy-metadraft]) or directly clone a repository (requires git to be installed).

Clone the cbmpy-metadraft git repository:

`git clone https://github.com/SystemsBioinformatics/cbmpy-metadraft.git`

MetaDraft ships with a single model template, to install the current template set, based on the BiGG2 repository, change into the MetaDraft directory and activate the template submodule:

`cd cbmpy-metadraft`
`git submodule update --init --remote -- modeldb/2018-1/`

You have now downloaded the latest version of MetaDraft, next please run the system test script `python systemtest.py` to check your installation and then proceed to install the required Python dependencies. How you do this will depend on whether you use Conda or PIP.

## Creating a custom Anaconda Python environment

The following commands will create custom virtual environments with all the Python dependencies needed to run MetaDraft. It is also possible to install the individual packages using "conda" or "pip" see the *requirements.txt* file for details of the required packages (see "pip" instructions below).

Adding the the core packages (minimal) to an existing CONDA environment:

`conda install -c bgoli -c sbmlteam sip pyqt python-libsbml xlrd xlwt cbmpy biopython`

### New Python 2 CONDA environment (recommended)

Ubuntu (Linux-64)

`conda env create -f conda_env_build_ubuntu.yml
conda activate metadraft2
./run.sh`

Windows 10 (Win-64)

`conda env create -f conda_env_build_win.yml
conda activate metadraft2
./runwin.bat`

### New Python 3 CONDA environment

Windows 10 (Win-64)

`conda env create -f conda_env_build_win_py3.yml
conda activate metadraft2
./runwin.bat`

## Using PyPI to install Python dependencies

### For Python 3

Using PyPI it is only possible to use Python 3. You can either install the packages useing the requirements.txt file (recommended):

`pip install -r requirements.txt`

or individually install the core packages (minimal):

`pip install sip PyQt5 python-libsbml xlrd xlwt cbmpy biopython`

Brett G. Olivier (b.g.olivier@vu.nl)

