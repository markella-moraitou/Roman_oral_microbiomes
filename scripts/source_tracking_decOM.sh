#!/bin/bash -l

#SBATCH -n 100
#SBATCH -p pelle
#SBATCH -t 5:00:00
#SBATCH -J decOM

# Run source tracking with decOM

#### SET ENVIRONMENT ####

[ ! -z $SLURM_SUBMIT_DIR ] && cd $SLURM_SUBMIT_DIR

echo "$SLURM_JOB_NAME"

DATADIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/results/samtools/filter 
OUTDIR=${PWD}/../decOM

conda activate decOM

echo $(module list)

echo $(conda list)

#### PREPARE INPUT ####

# Define path to sources, which sould be downloaded before using decOM
source=${decom_path}/decOM_sources
echo "Sources in $source"

mkdir -p $OUTDIR # Create output directory 

# Get sinks.txt file to use with -p_sinks as well as keys files to use with -p_keys, using the sample list

[ -f $OUTDIR/sinks.txt ] && rm -f $OUTDIR/sinks.txt && echo "Removing existing sinks.txt file"

ls $DATADIR/*unmapped.fastq.gz | while read input
do
  s=$(basename ${input%.unmapped.fastq.gz})
  if [[ ! -f $input ]]
  then
    echo "File $input doesn't exist. Skipping"
  else
    echo "Creating input files for $s"
    echo $s >> $OUTDIR/sinks.txt
    echo "${s} : ${PWD}/${input}" > $OUTDIR/${s}.fof
  fi
done

#### RUN SOURCETRACKING ####

cd $OUTDIR 

rm -rf decOM_output
decOM -p_sources $source -p_sinks sinks.txt -p_keys . -mem 30GB -t 20 -o decOM_output
