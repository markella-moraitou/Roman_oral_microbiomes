#!/bin/bash -l

#SBATCH -A naiss2024-22-84
#SBATCH -p memory
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 128
#SBATCH -t 24:00:00
#SBATCH -J clas_martina_parallel
#SBATCH -e err/%x-%j.error
#SBATCH -o err/%x-%j.out
#SBATCH --mail-user=chenyu.jin@su.se
#SBATCH --mail-type=FAIL
#SBATCH --mem=1760GB

# Script by Chenyu Jin, Center for Palaeogenetics, Swedish Museum of Natural History, Stockholm University

outdir="/cfs/klemming/projects/snic/snic2022-6-144/CHENYU/seq_classify_dardel/martina/clas"; export outdir
list_fq="/cfs/klemming/projects/snic/snic2022-6-144/CHENYU/seq_classify_dardel/martina/all_samples"

# Check if kelsey_files exists and is readable
if [[ ! -r $list_fq ]]; then
    echo "Error: $list_fq not found or not readable."
    exit 1
fi

# Function to run Kraken2
run_kraken2() {
    local db_path=$1
    local old_label=$4
    local label=$2
    local reads=$3
    local data_from=$(dirname "$reads")
    local name=$(basename "$reads" | sed 's/.fastq.gz//g' | sed "s/_${old_label}//g" | sed "s/_unclas//g" | sed "s/.merged//")

    # Create a unique temporary directory
    #local temp_db_dir=$(mktemp -d)

    # Copy or link the database to the temporary directory
    #ln -s "$db_path"/* "$temp_db_dir/"

    # Run Kraken2
    singularity exec --bind "$data_from,$outdir" kraken2.1.3.sif kraken2 \
        --db "$db_path" \
        --report "$outdir/${name}_${label}_report.txt" \
        --report-minimizer-data \
        --gzip-compressed \
        --threads $threads \
        --output "$outdir/${name}_${label}_output.txt" \
        --classified-out "$outdir/${name}_${label}_clas.fastq" \
        --unclassified-out "$outdir/${name}_${label}_unclas.fastq" \
        --memory-mapping "$reads"

    # Compress output files
    pigz "$outdir/${name}_${label}_clas.fastq" && \
    pigz "$outdir/${name}_${label}_unclas.fastq" && \
    pigz "$outdir/${name}_${label}_output.txt"

    # Clean up temporary directory
    rm -rf "$temp_db_dir"
}

export -f run_kraken2

threads=1; export threads

# Microbial GTDB
k2_db="/cfs/klemming/pdc/software/dardel/sw-uppmax/data/Kraken2_data/prebuilt/k2_gtdb_genome_reps_20241109";export k2_db
label1="1mic_GTDB";export label1

if [ ! -d "/tmp/$(basename $k2_db)" ]; then
            mkdir -p "/tmp/$(basename $k2_db)"
fi

cp -r $k2_db/*.k2d /tmp/$(basename $k2_db)

# Use xargs to parallelize the while loop with 2 tasks
xargs -P 125 -I {} bash -c 'run_kraken2 "/tmp/$(basename $k2_db)" "$label1" "{}" ""' < $list_fq
rm -r /tmp/$(basename $k2_db)

# Microbial Refseq
k2_db="/cfs/klemming/pdc/software/dardel/sw-uppmax/data/Kraken2_data/prebuilt/k2_pluspfp_20241228";export k2_db
label2="2mic_pluspfp";export label2

if [ ! -d "/tmp/$(basename $k2_db)" ]; then
            mkdir -p "/tmp/$(basename $k2_db)"
fi

cp -r $k2_db/*.k2d /tmp/$(basename $k2_db)

find "$outdir" -name "*_${label1}_unclas.fastq.gz" | xargs -P 125 -I {} bash -c 'run_kraken2 "/tmp/$(basename $k2_db)" "$label2" "{}" "$label1"'

rm -r /tmp/$(basename $k2_db)


#3plant
k2_db="/cfs/klemming/projects/snic/snic2022-6-144/KRAKEN2_DATABASES/PLANTS"; export k2_db
label3="3plant";export label3

if [ ! -d "/tmp/$(basename $k2_db)" ]; then
                    mkdir -p "/tmp/$(basename $k2_db)"
fi

cp -r $k2_db/*.k2d /tmp/$(basename $k2_db)

find "$outdir" -name "*_${label2}_unclas.fastq.gz" | xargs -P 125 -I {} bash -c 'run_kraken2 "/tmp/$(basename $k2_db)" "$label3" "{}" "$label2"'

rm -r /tmp/$(basename $k2_db)

#4Mammals
k2_db="/cfs/klemming/projects/snic/snic2022-6-144/KRAKEN2_DATABASES/MAMMALS"; export k2_db
label4="4mammal";export label4

if [ ! -d "/tmp/$(basename $k2_db)" ]; then
            mkdir -p "/tmp/$(basename $k2_db)"
fi

cp -r $k2_db/*.k2d /tmp/$(basename $k2_db)

find "$outdir" -name "*_${label3}_unclas.fastq.gz" | xargs -P 125 -I {} bash -c 'run_kraken2 "/tmp/$(basename $k2_db)" "$label4" "{}" "$label3"'

rm -r /tmp/$(basename $k2_db)

echo "Processing complete."
