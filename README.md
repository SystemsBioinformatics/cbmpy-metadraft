# cbmpy-metadraft
CBMPy Metadraft: a flexible and extensible genome-scale model reconstruction tool


To install the latest version of MetaDraft:
 - clone the cbmpy-metadraft git repository 
 - git clone https://github.com/SystemsBioinformatics/cbmpy-metadraft.git
 - MetaDraft ships with a single model template, to install the current template set, based on the BiGG2 repository, activate the template submodule
 - git submodule update --init --remote -- modeldb/2018-1/

You have now downloaded the latest version of MetaDraft, next please install the required Python dependencies, How you do this will depend on whether you use Conda or PIP. 
If you have Anaconda this is the recommended way to start, we will create a custom Conda evironment:

Ubuntu (Linux-64)
 - conda env create -f conda_env_build_ubuntu.yml
 - conda activate metadraft2
 - ./run.sh
 
Windows 10 (Win-64)
 - conda env create -f conda_env_build_win.yml
 - conda activate metadraft2
 - ./runwin.bat