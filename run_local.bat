:set MYPATH=%~dp0win64\ncbi-blast-2.2.26
set MYPATH=%~dp0win64\ncbi-blast-2.2.26;%~dp0tools\win64\jre1.8.0_192\bin;%~dp0tools\win64\Perl64\bin
set "PATH=%PATH%;%MYPATH%"
echo PATH
python metadraft.py
