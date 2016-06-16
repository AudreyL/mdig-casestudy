#!/bin/bash

###############################################################################
# run mdig-notIntegratedInGrass.py from bash
###############################################################################
# Audrey Lustig
# Bio-Protection Research Centre, Lincoln University, NZ
# March 2016

homeDirectory='/home/audrey/Desktop/mdig-casestudy'
DataDirectory=$homeDirectory/L.Dispar/Data/SurvivalLayer/*/

for dir in $DataDirectory;
do
	cd $dir 
	num_rep=$(ls | wc -l)

	name="${dir%/*}"
	name_dir="${name#*Data/SurvivalLayer/}"

	ListAllFile='*'
	AsciiFiles=$dir$ListAllFile
	
	for f in $AsciiFiles;
	do

	filemane_asc="${f##*/}"
	
	python $homeDirectory/L.Dispar/mdig-notIntegratedInGrass.py $name_dir $filemane_asc $homeDirectory

	done
done
