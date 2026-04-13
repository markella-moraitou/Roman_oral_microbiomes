#!/bin/bash -l
#SBATCH -A snic2022-5-255
#SBATCH -p core
#SBATCH -n 8 
#SBATCH --array=2,3,11,12%4
#SBATCH -t 04:00:00
#SBATCH -J demux
#SBATCH --mail-type=FAIL

module load bioinfo-tools AdapterRemoval

# date
ymd=$(date +%y%m%d)
echo "$ymd"

threads=$SLURM_NTASKS

READDIR=${PWD}/../raw_reads/delivery06069/VD-3269/220617_A00605_0444_AHMGCKDSX3/

IN5=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${READDIR}/index_combs.txt | awk '{print $1}')
echo "$IN5"
IN7=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${READDIR}/index_combs.txt | awk '{print $2}')
echo "$IN7"
lane=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${READDIR}/index_combs.txt | awk '{print $3}')
echo "$lane"


OUTDIR=${PWD}/../demultiplex/delivery06069/0_demux_${ymd}/S$IN5/L00$lane
mkdir -p "$OUTDIR"

BCLIST=$OUTDIR/${IN5}_bclist.txt
grep "$IN5,$IN7" ${PWD}/../demultiplex/delivery06069/0_demux_220705/human_barcode_demux_280622.txt | tr "," "\t" | awk '{print $1 "\t" $2 "\t" $3}' > $BCLIST

BC5DIR=${READDIR}/Sample_VD-3269-Index${IN5}
BC7DIR=${READDIR}/Sample_VD-3269-Index${IN7}

# the string {name} in the name of the output file and give each adapter a name.
FWD=$(find "${BC5DIR}"/VD-3269-Index${IN5}_S13${IN5}_L00${lane}_R1_001.fastq.gz)
REV=$(find "${BC7DIR}"/VD-3269-Index${IN7}_S13${IN7}_L00${lane}_R2_001.fastq.gz)

# trimns: If set, trim ambiguous bases (N) at 5'/3' termini [default: off]
cd $OUTDIR
AdapterRemoval \
  --file1 $FWD \
  --file2 $REV \
  --demultiplex-only \
  --barcode-mm 2 \
  --gzip \
  --basename "$OUTDIR"/ \
  --threads $threads \
  --barcode-list "$BCLIST" \
  --mate-separator ":"