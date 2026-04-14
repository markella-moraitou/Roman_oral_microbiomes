#!/bin/bash

#SBATCH -p core
#SBATCH --ntasks 8 --mem-per-cpu=1G
#SBATCH -t 02:00:00
#SBATCH -J adapterremovalcount

# USAGE: sbatch count_adapterremoval.sh [DIRNAME]
# INPUT: [DIRNAME], contaning .pe.settings files with AdapterRemoval info
# OUTPUT: Number_fulllength_collapsed.txt in your [DIRNAME]

module load gnuparallel

ymd=$(date +%y%m%d)
echo $ymd

threads=8

find $1/*.pe.settings | parallel -j $threads "echo {}; grep "Number of full-length collapsed pairs" {}" >> $1/Number_fulllength_collapsed.txt