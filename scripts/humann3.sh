#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 5-00:00:00
#SBATCH -J HUMAnN3

# Humann3 pipeline with Metaphlan2 for tax profiling and mapping/DIAMOND for functional assignments

#### SET ENVIRONMENT ####

echo "$SLURM_JOB_NAME"

echo "Number of tasks: ${SLURM_NTASKS}"

DATADIR=${PWD}/../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/results/samtools/filter 
OUTDIR=${PWD}/../HUMANn3
CONDAPATH=$(conda info --base)/envs/humann_env
MAPPING=${CONDAPATH}/databases/utility_mapping/

#conda create -n humann_env -c bioconda -c conda-forge python=3.9 metaphlan=3.0 bowtie2 diamond pip
conda activate humann_env
#pip install humann==3.9
#metaphlan --install --index mpa_v30_CHOCOPhlAn_201901

#Install database
#cd $CONDAPATH
#mkdir -p databases
#cd databases

#humann_databases --download chocophlan full .
#humann_databases --download uniref uniref90_diamond .
#humann_databases --download utility_mapping full .

echo $(module list)

echo $(conda list)

#### RUN HUMANn3 ####

cd $DATADIR || exit
ls *unmapped.fastq.gz | while read i 
do
    if [[ -f $OUTDIR/${i%.unmapped.fastq.gz}_genefamilies.tsv ]]
    then
        continue
    else
      humann \
      --nucleotide-database ${CONDAPATH}/databases/chocophlan/ \
      --protein-database ${CONDAPATH}/databases/uniref/ \
      --metaphlan ${CONDAPATH}/databases/metaphlan \
      --metaphlan-options "--bowtie2db ${CONDAPATH}/lib/python3.9/site-packages/metaphlan/metaphlan_databases/" \
      --input "$i" \
      --output $OUTDIR \
      --output-basename "${i%.unmapped.fastq.gz}" \
      --threads "$SLURM_NTASKS"
    fi
done

#regroup gene families into different functional categories (KEGG orthogroups and GO terms)
#renormalization
#joining tables

cd $OUTDIR

for i in *genefamilies.tsv;
do
    #Get KO groups
    humann_regroup_table --input $i --output ${i%genefamilies.tsv}KOgroups.tsv --custom $MAPPING/map_ko_uniref90.txt.gz
    humann_renorm_table --input ${i%genefamilies.tsv}KOgroups.tsv --output ${i%genefamilies.tsv}KOgroups_cpm.tsv --units cpm
    #Normalize path abundances
    #humann_renorm_table --input ${i%genefamilies.tsv}pathabundance.tsv --output ${i%genefamilies.tsv}pathabundance_cpm.tsv --units cpm;
    #Get GO terms
    humann_regroup_table --input $i --output ${i%genefamilies.tsv}go.tsv --custom $MAPPING/map_go_uniref90.txt.gz
    humann_renorm_table --input ${i%genefamilies.tsv}go.tsv --output ${i%genefamilies.tsv}go_cpm.tsv --units cpm
done

# join all normalized output files into one table
humann_join_tables --input $OUTDIR --output all_KOgroups_cpm.tsv --file_name KOgroups_cpm
humann_join_tables --input $OUTDIR --output all_GOterms_cpm.tsv --file_name go_cpm
#humann_join_tables --input $OUTDIR --output all_pathabund_cpm.tsv --file_name pathabundance_cpm
#humann_join_tables --input $OUTDIR --output all_pathcov.tsv --file_name pathcoverage

# get unstratified info
#grep "|" -v $OUTDIR/all_KOgroups_cpm.tsv > $OUTDIR/all_KOgroups_cpm_unstratified.tsv
#grep "|" -v $OUTDIR/all_GOterms_cpm.tsv > $OUTDIR/all_GOterms_cpm_unstratified.tsv
#grep "|" -v $OUTDIR/all_pathabund_cpm.tsv > $OUTDIR/all_pathabund_cpm_unstratified.tsv
#grep "|" -v $OUTDIR/all_pathcov.tsv > $OUTDIR/all_pathcov_unstratified.tsv


#get human_readable names
#humann_rename_table --input all_GOterms_cpm_unstratified.tsv -n go --output all_GOterms_cpm_unstratified_renamed.tsv
humann_rename_table --input all_GOterms_cpm.tsv -n go --output all_GOterms_cpm_renamed.tsv
