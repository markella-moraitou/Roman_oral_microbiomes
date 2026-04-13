#!/bin/bash
#SBATCH -A snic2021-5-477
#SBATCH -p core
#SBATCH -n 4
#SBATCH --mem-per-cpu 2GB
#SBATCH -t 0:30:00
#SBATCH --array=1-96%15
#SBATCH -J bracken

# date
ymd=$(date +%y%m%d)

# modules
module load bioinfo-tools Kraken Kraken2 Bracken
echo $(module list)
echo $(readlink -f $KRAKEN2_DEFAULT_DB)

# Bracken Settings
kmer=35
read_length=65
threads=$SLURM_CPUS_ON_NODE

BASEDIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022
KRAKENDIR=$BASEDIR/results/metagenomic_classification/kraken/standard

# Names from the sample list
i=`find $KRAKENDIR/*.unmapped.fastq.kraken2_report | sed -n "$SLURM_ARRAY_TASK_ID"p |  awk '{print $1}'`
name=`echo $i |  sed 's!.*/!!' | cut -d . -f 1` # file name

OUTDIR=$BASEDIR/results/metagenomic_classification/bracken/standard
mkdir -p $OUTDIR

MY_DB=/proj/sllstore2017021/nobackup/SHARED_REFERENCES/KRAKEN2/standard

# bracken
# -d MY_DB -i INPUT -o OUTPUT -r READ_LEN -l LEVEL -t THRESHOLD
# no threshold i.e. 0 reads required before reestimation (-t 0)
# will apply abundance threshold independently

cd $OUTDIR || exit
bracken -i $i -o $OUTDIR/${name}_kraken2_report_bracken.txt -d $MY_DB -r $read_length -l S -t 0
