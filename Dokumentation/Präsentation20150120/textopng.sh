#!/bin/bash
read -p "Dateinummer?:" nummer
Datei=folie$nummer
echo $Datei
latex $Datei.tex
dvipng -T tight -D 512 $Datei.dvi
rm *.log *.eps *.aux *.dvi
