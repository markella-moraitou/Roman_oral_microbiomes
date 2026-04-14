#!/bin/bash -l

#SBATCH -n 20
#SBATCH -p pelle
#SBATCH -t 3:00:00
#SBATCH -J mapdamage

# Mapping metagenomics reads to microbial genomes and running MapDamage

#### SET ENVIRONMENT ####

[ ! -z $SLURM_SUBMIT_DIR ] && cd $SLURM_SUBMIT_DIR

conda activate base
module load pigz BWA mapDamage

echo $SLURM_JOB_NAME
echo $(module list)

if [ -z "$SLURM_NTASKS" ] ; then
  processes=1
else
  processes=$SLURM_NTASKS
fi

BASEDIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022
DATADIR=$BASEDIR/results/metagenomic_classification/kraken/standard
SEQDIR=$BASEDIR/results/adapterremoval/output/
OUTDIR=${PWD}/../mapdamage

mkdir -p $OUTDIR

# Download references for the taxa of interest
mkdir -p $OUTDIR/ref_genomes

[ ! -f $OUTDIR/ref_genomes/genome_summary.csv ] && echo "OTU,species,genome_file,genome_type" > $OUTDIR/ref_genomes/genome_summary.csv

tail -n+2 ../output/mapdamage_taxa_samples.csv | while IFS=, read -r taxid species _ _
do
  if [ ! -f $(find . -name "*.fna" | head -1 ) ]
  then
      cd $OUTDIR/ref_genomes
      new_species=$(echo $species | tr " " _)
      echo "Downloading reference genome for $species"
      datasets download genome taxon $taxid --reference --filename ${new_species}_${taxid}.zip
      type="reference"
      # If no reference genomes was found
      if [ ! -f ${new_species}_${taxid}.zip ]
      then
          echo "Reference genome for $species not found, downloading complete assembly"
          datasets download genome taxon $taxid --filename ${new_species}_${taxid}.zip
          type="not_reference"
      fi
      unzip ${new_species}_${taxid}.zip 
      genome=$(find ncbi_dataset -name "*.fna" | head -1 ) # Pick one genome
      mv $genome ${new_species}_${taxid}_$(basename ${genome})
      rm -rf ncbi_dataset ${new_species}_${taxid}.zip md5sum.txt README.md
      echo "$taxid,${species},$(basename ${genome}),${type}" >> genome_summary.csv
      cd -
  fi
done

echo taxid,species,sample,mapped_reads > $OUTDIR/mapped_reads.csv

# Map reads to the taxon genome from the sample with the highest abundance
tail -n+2 ../output/mapdamage_taxa_samples.csv | while IFS=, read -r taxid species sample _
do  
    cd $OUTDIR
    new_species=$(echo $species | tr " " _)
    # Identify Kraken and fastq inputs
    krakenin=${DATADIR}/${sample}.unmapped.fastq.kraken.out
    # Merge two lanes of adapterremoval output 
    cat ${SEQDIR}/${sample}_L00*_F.fastq.pG.fq_L*.pe.combined.fq.gz > ${sample}_adapterremoval.fa.gz
    seqin=${sample}_adapterremoval.fa.gz
    ## Index and map
    ref=$(find $OUTDIR/ref_genomes -name ${new_species}_${taxid}*.fna)
    alignment=$OUTDIR/${sample}_${new_species}_alignment.sam
    if [ ! -f $ref ]
    then
        echo "Genome not found for $species. Skipping..."
        continue
    fi
    #bwa index -a bwtsw ${ref}
    bwa mem -t $processes ${ref} ${seqin} > $alignment
    samtools view $alignment -F 4 -@ 10 -O bam > ${alignment%.sam}.bam
    nreads=$(samtools view -c ${alignment%.sam}.bam )
    echo "$taxid,$species,$sample,$nreads" >> $OUTDIR/mapped_reads.csv
    ## Run mapdamage
    mapDamage -i ${alignment%.sam}.bam -r ${ref} -d mapdamage_${sample}_${new_species}
    cd -
done
