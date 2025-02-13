# InParanoid Diamond update

In 2022 InParanoid was reimplemented using DIAMOND and released under a much more usable GPLv3 licence. Finally Metadraft can be reimplemented as a complete OS package.

- InParanoid-DIAMOND Publication: https://academic.oup.com/bioinformatics/article/38/10/2918/6561543
- InParanoid-DIAMOND code: https://bitbucket.org/sonnhammergroup/inparanoid/src/master/

## InParanoid-DIAMOND  - Ubuntu 24.04.2 LTS
*Mostly extracts from the InParanoid-DIAMOND README file.*

First the dependencies - please check that you know what these are before installing

```bash
sudo apt install bioperl libmoose-perl libparallel-forkmanager-perl
````

Next grab the inParanoid-DIAMOND source

```bash
git clone https://bitbucket.org/sonnhammergroup/inparanoid.git
```

## Installing DIAMOND

```bash
mkdir diamond
cd diamond
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```

## Testing inParanoid-DIAMOND
If everything has gone according to plan this test should work, go to the inparanoid directory which should be next to the diamond directory and execute:

```bash
perl inparanoid.pl -diamond-path ../diamond/diamond -input-dir ./testInput/

Starting Diamond searches...
DIAMOND searches took 0.05 seconds
Done with DIAMOND searches. Starting ortholog detection...
Reading and sorting homologs took 0.01 seconds
Finding and sorting orthologs took 0.00 seconds
Reading paralogous hits took 0.00 seconds
Finding in-paralogs took 0.00 seconds
mysql output saved to ./output/SQLtable.EC:pid30261-SC:pid30261
Finding bootstrap values and printing took 0.00 seconds
Running InParanoid for EC:pid30261 - SC:pid30261 took 1.00 seconds

Total real time: 1.00 seconds
Total user time: 2.82 seconds
Total system time: 0.14 seconds
```

If you get a DIAMOND not in path error then adjust the -diamond-path flag for your system.
