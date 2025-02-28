#!/bin/sh

#set JRE_PATH=c:\DevTools\jdk-23.0.2\bin
#set PERL_PATH=c:\DevTools\strawberry\perl\bin;c:\DevTools\strawberry\c\bin
#set DIAMOND_PATH=%~dp0diamond
##set INPARANOID_PATH=%~dp0inparanoid
#set "PATH=%PATH%;%JRE_PATH%;%PERL_PATH%;%DIAMOND_PATH%;%INPARANOID_PATH%"

export CDIR=$(pwd)
export DIAMOND_PATH=$CDIR/diamond
export INPARANOID_PATH=$CDIR/inparanoid
export PATH=$PATH:$CDIR:$DIAMOND_PATH

echo " "
echo # Metadraft system path
echo $PATH
echo " "

python metadraft.py

