@echo off

where formatdb.exe /q
if %errorlevel% EQU 0 (
  echo BLAST OK
) else (
  set path="%path%";"%~dp0\tools\win64\ncbi-blast-2.2.26\"
  ::echo "Work path %~dp0"
  echo Using builtin BLAST "%~dp0\tools\win64\ncbi-blast-2.2.26\"
)

python metadraft.py

