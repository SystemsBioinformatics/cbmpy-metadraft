#!/bin/bash
cd /inparanoid/
export PATH=/inparanoid/:$PATH

if [[ "$@" == "-h" || "$@" == "-help" ]]
then
perl inparanoid.pl -h
elif [[ "$@" == *"-f1"* || "$@" == *"-f2"* || "$@" == *"-input-dir"* ]]
then
echo -e "\nPlease dont use the -f1, -f2 or -input-dir parameters when running the container version of InParanoid. When running the docker container you have to mount your input dir to the container, and leave the input-file-parameters (-f1, -f2 and -input-dir) empty. \n\nUse the following command to run the container version: \ndocker run -v <path/to/your/input/files>:/input/ -v <path/where/you/want/output/files>:/output/ sonnhammer/inparanoid\n"
elif [[ "$@" == *"-out-dir"* ]]
then
echo -e "\nPlease dont use the -out-dir parameter when running the container version of InParanoid. When running the docker container you have to mount your output dir to the container, and leave the out-dir parameter empty. \n\nUse the following command to run the container version: \ndocker run -v <path/to/your/input/files>:/input/ -v <path/where/you/want/output/files>:/output/ sonnhammer/inparanoid\n"
else
cd /output
cp /inparanoid/inparanoid.pl .
cp /inparanoid/diamondParser.pl .
cp /inparanoid/blast_parser.pl .
cp /inparanoid/seqstat.jar .
cp /inparanoid/help .
cp -rp /inparanoid/matrices .
perl inparanoid.pl -input-dir /input/ -out-dir /output/ "$@"
rm ./inparanoid.pl
rm ./diamondParser.pl
rm ./blast_parser.pl
rm ./seqstat.jar
rm ./help
rm -R ./matrices
fi





