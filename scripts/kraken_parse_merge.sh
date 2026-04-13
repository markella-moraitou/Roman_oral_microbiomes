#!/bin/bash

module load python3

BASEDIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022
KRAKENDIR=$BASEDIR/results/metagenomic_classification/kraken/standard
BRACKENDIR=$BASEDIR/results/metagenomic_classification/bracken/standard

krakenfiles=$(find $KRAKENDIR/*.kraken2_report)

# using Adrian's version of parse and merge (which as extra parameters for bracken formatted output)
# kraken output

# first remove old csv files
rm -f $KRAKENDIR/*.csv

for i in $krakenfiles
do  
  name=$(basename $i .kraken2_report)
  python3 kraken_parse.py -c 1 -or $KRAKENDIR/${name}.read_kraken_parsed.csv -ok $KRAKENDIR/${name}.kmer_kraken_parsed.csv $i
done

cd ${KRAKENDIR} || exit

python3 merge_kraken_res.py -or read_kraken_parsed.csv -ok kmer_kraken_parsed.csv

# bracken output
# remeber to not include output for kmer files, not possible with bracken output
brackenfiles=$(find $BRACKENDIR/*_bracken_species.kraken2_report)
# first remove old csv files
rm -f $BRACKENDIR/*.csv

for i in $brackenfiles
do  
  name=$(basename $i .kraken2_report)
  python3 kraken_parse.py -c 1 -or $BRACKENDIR/${name}.read_bracken_parsed.csv $i --bracken
done

cd ${BRACKENDIR} || exit
python3 merge_kraken_res.py -or read_bracken_parsed.csv --bracken
