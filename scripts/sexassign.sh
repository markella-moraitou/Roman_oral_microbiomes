#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J sexassign

# Determine genetic sex using sexassign https://github.com/grahamgower/sexassign
# Gower et al. (2019), https://doi.org/10.1073/pnas.1903275116

module load bioinfo-tools samtools bbmap

echo $SLURM_JOB_NAME
echo $(module list)

OUTDIR=${PWD}/../sexassign
DATADIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/results/samtools/filter/
SCRIPTDIR=${PWD}

cd $DATADIR

# Chromosomes to keep (ignore Y chr and all unplaced scaffolds)

chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

# Initiate file to record number of reads
echo "sample,samtools_filtered,chr_filtered" > ${OUTDIR}/readcount_sexassign.csv

# Use BAM files from mapping to the human genome
# downsample
for i in *.bam
do
  sample=${i%_PE.mapped_chr11.bam}
  # Extract only full chromosomes. Need to do this one by one and then merge them
  for chr in ${chromosomes[@]}
  do
     samtools view -b -o ${OUTDIR}/${i%.bam}_${chr}.bam $i $chr
  done
  samtools merge -b <( ls ${OUTDIR}/${i%.bam}_chr*.bam ) -o ${OUTDIR}/${i%.bam}_filter_chr.bam
  rm -f ${OUTDIR}/${i%.bam}_chr*.bam
  
  pre_rc=$(samtools view -c $i)
  post_rc=$(samtools view -c ${OUTDIR}/${i%.bam}_filter_chr.bam)
  echo "$sample,$pre_rc,$post_rc"  >> ${OUTDIR}/readcount_sexassign.csv
done

cd $OUTDIR

# Get coverage statistics -- need to index first
for i in *filter_chr.bam
do
  samtools index -b $i
  samtools idxstats ${i} > ${OUTDIR}/${i%.bam}.idxstats
done

# Run sexassign with different minimum read count thresholds

python ${SCRIPTDIR}/sexassign.py *.idxstats -w -o sexes.pdf > sexes.txt

python ${SCRIPTDIR}/sexassign.py *.idxstats -w -o sexes_1000.pdf -n 1000 > sexes_1000.txt

python ${SCRIPTDIR}/sexassign.py *.idxstats -w -o sexes_100.pdf -n 100 > sexes_100.txt
