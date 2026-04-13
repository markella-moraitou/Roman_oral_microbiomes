#!/bin/bash

#SBATCH -A snic2021-5-477
#SBATCH -p core
#SBATCH --ntasks 8 --mem-per-cpu=1G
#SBATCH -t 02:00:00
#SBATCH -J fqcount
#SBATCH -o /proj/sllstore2017021/nobackup/MARTINA/logs/slurm-%A_%a.out

# USAGE: sbatch quick_count.sh [DIRNAME]
# INPUT: [DIRNAME], contaning fastq files
# OUTPUT: fastq_rc.txt in your [DIRNAME]

module load gnuparallel

ymd=$(date +%y%m%d)
echo $ymd

threads=8

find $1/*.fq.gz | parallel -j $threads "echo {}; zgrep -c "^@" {}" >> $1/fastq_rc.txt

