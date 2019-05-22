# cbmpy-metadraft
CBMPy Metadraft: a flexible and extensible, GUI-based genome-scale model reconstruction tool that supports multiple Systems Biology standards.

## **Important** upgrading older versions of MetaDraft to version 0.9.2

MetaDraft version 0.9.2 now defaults to the *2019-1* template library, if you are working directly in a cloned repository please initialise the new template library:

```Shell
cd cbmpy-metadraft
git submodule update --init --remote -- modeldb/2019-1/
```

If you want to switch back to the *2018-1* library and results please edit the `_metadraft.cfg` file and change the value  of the `"metadraft_db_version"` key to `"2018-1"` (it is currently set to `"2019-1"`).


## Installation and system dependencies

NEW: The MetaDraft User Guide is now available online and for download from the MetaDraft [website](https://systemsbioinformatics.github.io/cbmpy-metadraft/).

### System dependencies

MetaDraft requires working versions of Python, Perl, Java and NCBI Blast2. Once you have downloaded MetaDraft and it's Python dependencies, please run the *system test* script to that can be found in the MetaDraft directory `python systemtest.py`.

Typically the default Java runtime (JRE) installed on your system is all that is require, for Windows this can be obtained from [https://www.java.com/]. On Linux please use your favourite package manager, for example, when using Ubuntu/Debian based oeprating systems, try `sudo apt-get install default-jre`. 

MetaDraft has been successfully tested using Windows 10, Ubuntu Linux 16.04/18.04, Python 2.7, Python 3.6, PyQt4 and PyQt5. MetaDraft is developed on Windows Anaconda Python 2.7 with PyQt5.

### Getting MetaDraft

To install the latest version of MetaDraft, download from GitHub ([https://github.com/SystemsBioinformatics/cbmpy-metadraft]) or directly clone a repository (requires git to be installed).

Clone the cbmpy-metadraft git repository:

`git clone https://github.com/SystemsBioinformatics/cbmpy-metadraft.git`

### Getting and installing the MetaDraft template library

MetaDraft ships with a single model template, if you have cloned the repository you can activate the current template set, based on the BiGG2 repository, change into the MetaDraft directory and activate the template submodule:

`cd cbmpy-metadraft`

`git submodule update --init --remote -- modeldb/2019-1/`

If you have downloaded master repository as a zip file or as a release, download the `template-models-2018-1.zip` template archive provided with the [latest release](https://github.com/SystemsBioinformatics/cbmpy-metadraft/releases). Unzip this archive into the `metadraft/modeldb` directory.

MetaDraft is archived on Zenodo, please see CITATION.cff for more details [![DOI](https://zenodo.org/badge/132483758.svg)](https://zenodo.org/badge/latestdoi/132483758)

## Creating a custom Anaconda Python environment

You have now downloaded the latest version of MetaDraft, next please run the system test script `python systemtest.py` to check your installation and then proceed to install the required Python dependencies. How you do this will depend on whether you use Conda or PIP.

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

## Software citation

Please cite this software if you use it (the author(s) enthusiasm to carry on maintaining this project depend on it), see CITATION.cff for more details.

```text
 Brett G. Olivier. (2018, December 18). 
 SystemsBioinformatics/cbmpy-metadraft: Metadraft. 
 Zenodo. http://doi.org/10.5281/zenodo.2398336
```




Brett G. Olivier (b.g.olivier@vu.nl)

