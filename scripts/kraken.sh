#!/bin/bash
#SBATCH -A snic2021-5-477
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -J kraken
#SBATCH -M snowy
#SBATCH --array=1-20%20
#SBATCH -o /proj/sllstore2017021/nobackup/MARTINA/logs/slurm-%A_%a.out
#SBATCH --mail-type=FAIL

# Microbial assignments using kraken
# genbank nt database
# build by UPPMAX $KRAKEN_DB (updated started of each month)

module load bioinfo-tools Kraken2
echo $(module list)
echo $(readlink -f $KRAKEN2_DEFAULT_DB)

#### SET VARIABLES
ymd=$(date +%y%m%d)
echo $ymd

threads=$SLURM_CPUS_ON_NODE
echo $threads

BASEDIR=${PWD}../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022
FASTQDIR=$BASEDIR/results/samtools/filter

#### KRAKEN
## Set up database
MY_DB_DIR=/proj/sllstore2017021/nobackup/SHARED_REFERENCES/KRAKEN2/standard

#### FIND FILES, GENERATE NAME VARIABLES
# Names from the sample list
# remove reindeer from list, already analysed with the current db
i=$(find $FASTQDIR/*.unmapped.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p)
name=$(echo $i | sed 's!.*/!!' | cut -d . -f 1) # file name

# set dirs
OUTDIR=${BASEDIR}/kraken/standard
mkdir -p $OUTDIR

kraken2 --db $MY_DB_DIR $i --threads $threads --report-minimizer-data --report $OUTDIR/${name}.unmapped.fastq.kraken2_report --output $OUTDIR/${name}.unmapped.fastq.kraken.out