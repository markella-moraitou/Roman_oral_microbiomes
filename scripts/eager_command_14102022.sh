# export SNIC_TMP="/scratch/martinaf"
# module load bioinfo-tools nf-core nf-core-pipelines/latest Nextflow/22.10.0
# export NXF_HOME="../nf-core-eager"

nextflow run nf-core/eager -r 2.4.4 -profile uppmax -work-dir ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/work \
        --input ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/RDC_samplemap_14102022.tsv \
        --fasta /sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
        --bwa_index /sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/ \
        --large_ref true --save_reference true \
        --outdir ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/results --complexity_filter_poly_g true --clip_readlength 35 --clip_min_read_quality 30 --mergedonly true --run_post_ar_trimming true --bt2n 0 --run_bam_filtering true --bam_mapping_quality_threshold 30 --bam_filter_minreadlength 35 --bam_unmapped_type fastq --run_trim_bam true --run_genotyping true --genotyping_tool angsd --genotyping_source trimmed --run_vcf2genome false --run_multivcfanalyzer false --write_allele_frequencies true --run_mtnucratio true --mtnucratio_header chrM --run_sexdeterrmine false --run_nuclear_contamination true --contamination_chrom_name chrX --run_metagenomic_screening true --metagenomic_complexity_filter false --metagenomic_tool kraken --database /proj/sllstore2017021/nobackup/SHARED_REFERENCES/KRAKEN2/standard/ 

#to get genotyping output 17/02/2023
nextflow run nf-core/eager -r 2.4.4 -profile uppmax -work-dir ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/work \
         --input ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/RDC_samplemap_14102022.tsv \
         --fasta /sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
         --bwa_index /sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/ \
         --large_ref true --save_reference true \
         --outdir ../nf-core-eager/human_wholedataset_mergfalse_nosexdet_14102022/results \
         --skip_adapterremoval true \
         --complexity_filter_poly_g true --clip_readlength 35 --clip_min_read_quality 30 --mergedonly true \
         --run_post_ar_trimming true \
         --bt2n 0 \
         --run_bam_filtering true --bam_mapping_quality_threshold 30 --bam_filter_minreadlength 35 --bam_unmapped_type fastq --run_trim_bam false \
         --run_genotyping true --genotyping_tool ug --genotyping_source raw --gatk_call_conf 30 --gatk_ploidy 2 --gatk_downsample 250 \
         --run_vcf2genome true  --run_multivcfanalyzer true --write_allele_frequencies true \
         --run_mtnucratio false --run_sexdeterrmine false --run_nuclear_contamination false --run_metagenomic_screening false