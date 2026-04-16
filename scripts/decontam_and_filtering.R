library(tidyverse)
library(vegan)
library(data.table)
library(phyloseq)
library(FEAST)
library(decontam)
library(microbiome)
library(cuperdec)
library(decontam)
library(microbiomeutilities)
library(cowplot)
library(scales)

# SET THEME AND PALETTES ------------------------------------

sample_palette <- c("Human" = "#3DA417", "Soil" = "#753200", "Ext Blank" = "#C6C3C1", "LP Blank" = "#F3F2F1")

necr_palette <- c("Lucus Feroniae" = "#0D1F61", "Selvicciola" = "#DA8434", "Isola Sacra" = "#2AAE2A")

health_palette <- c("Healthy" = "#5EABAF", "Healthy in caries mouth" = "#DAA334", "Caries" = "#B51C18")

source_palette <- c("oral"="#CC121B", caries = "#7F0006", "gut"="#A97004", "soil"="#3DC02F", "skin"="#EBAD39", "labcontam"="#32569C", "unknown"="#CBCBCB")

source_palette2 <- c("aOral"="#CC121B", "mOral" = "#D46A6A", "Sediment.Soil"="#3DC02F", "Skin"="#EBAD39", "Unknown"="#CBCBCB")


wd <- getwd()

# IMPORT DATA -------------------------------------------------

# SAMPLE DATA
meta <- readxl::read_xlsx(file.path("..", "R_input", "Masterlist_RDC.xlsx"), sheet = 1) %>% filter(!is.na(LP_ID)) %>% as.data.frame()
rownames(meta) <- meta$LP_ID

# Get better names
colnames(meta) <- make.names(colnames(meta)) %>% gsub(pattern="..", replacement=".", fixed = TRUE)

# Get some important variables as factors so that they get ordered correctly in plots
meta <- meta %>%
  mutate(Sample.Category = factor(Sample.Category, levels = c("Human", "Soil", "Ext Blank", "LP Blank"))) %>%
  mutate(Necropolis = factor(Necropolis, levels = c("Isola Sacra", "Lucus Feroniae", "Selvicciola"))) %>%
  mutate(Tooth.health = factor(Tooth.health, levels = c("Healthy", "Healthy in caries mouth", "Caries")))

# OTU TABLE
k.mcat <- read_csv(file.path("..", "R_input", "read_bracken_parsed_human.csv")) %>%
  column_to_rownames("TAXID")

# SOURCE OTU TABLE FOR FEAST
otu_sources <- read.csv(file.path("..", "R_input", "read_bracken_parsed_xfeast.csv"), header = T,comment.char = "",sep=",") %>%
  filter(TAXID %in% rownames(k.mcat)) # Keep only taxa that exist in sink OTU table

# decOM results
decom <- read.csv(file.path("..", "R_input", "decOM_output.csv"), header = T,comment.char = "",sep=",")

# Lists of oral taxa
# HOMD
homd.df <- read_delim(file.path("..", "R_input", "HOMD_taxon_table2022-08-12_1660320566.txt"), col_names = T, delim = "\t",skip=1) %>%
  dplyr::rename(TAXID=NCBI_taxon_id) %>% 
  distinct(TAXID, .keep_all = TRUE) %>%
  # Keep only oral taxa
  filter(grepl("Oral$", Body_site))

# James's list of hominid core oral microbiome taxa.
hominid.df <- read.csv(file.path("..", "R_input", "fellowyates2021-taxa-wID.csv"), header = T) %>%
  dplyr::rename(TAXID=tax_id)

# common contaminants
# only genus level available, so going to have to match by name
contaminant_list <- read.csv(file.path("..", "R_input", "common-lab-contaminants.csv")) %>%
  dplyr::rename(TAXID=tax_id)

# GET TAXONOMIC INFO -------------------------------------------------

# TAXONOMIC DATA
if (!file.exists(file.path("..", "R_output", "microbial_taxonomy.csv"))) {
  full_taxonomy <- list()
  i = 0
  taxids <- rownames(k.mcat)
  for (id in taxids) {
    i = i + 1
    tax <- tryCatch({taxize::classification(id, db="ncbi")},
                         # Try to catch the errors since this command often gives HTTP 400 errors
                         error = function(e) {cat("Error occurred for taxon ", id, ":")
                           cat(conditionMessage(e))
                           cat("...Retrying...")
                           Sys.sleep(1)
                           taxize::classification(id, db="ncbi")
                         })
    # Print info
    inf <- tax[[1]][nrow(tax[[1]]), 1]
    cat(i, "ID ", id, ": ", inf, "\n")
    full_taxonomy[[id]] <- tax[[1]]
  }
  
  # Get taxonomic info in a table format
  extract_taxonomy_info <- function(taxonomy) {
    desired_ranks <- c("id", "species", "genus", "family", "order", "class", "phylum")
    
    id <- taxonomy$id[nrow(taxonomy)]
    
    # Create a table with ranks
    tax_info <- bind_rows(taxonomy, data.frame(name = id, rank = "id")) %>%
      filter(rank %in% desired_ranks) %>%
      select(rank, name) %>%
      pivot_wider(names_from = "rank", values_from = "name")
    
    return(tax_info)
  }
  
  # Get taxonomic info in a column and reorganise
  taxonomy_mat <- as.matrix(bind_rows(lapply(full_taxonomy, extract_taxonomy_info)))
  rownames(taxonomy_mat) <- taxonomy_mat[,"id"]
  taxonomy_mat <- taxonomy_mat[, -which(colnames(taxonomy_mat) == "id")]
  # Save output
  saveRDS(full_taxonomy, file = file.path("..", "R_output", "microbial_taxonomy_full.RDS"))
  write.csv(taxonomy_mat, file = file.path("..", "R_output", "microbial_taxonomy.csv"), quote = FALSE, row.names = TRUE)
} else {
  taxonomy_mat <- read.csv(file = file.path("..", "R_output", "microbial_taxonomy.csv")) %>% column_to_rownames("X") %>% as.matrix
}

# CREATE PHYLOSEQ OBJECTS -------------------------------------------------

# Combine into a phyloseq object
k.ps <- phyloseq(otu_table(k.mcat, taxa_are_rows = TRUE), sample_data(meta), tax_table(taxonomy_mat))

# update metadata with bracken counts/richness
k.ps@sam_data$Bracken.readct.total <-  stack(sample_sums(k.ps))[,"values"]

# Only keep species level assignments
k.ps.species <- k.ps %>% subset_taxa(!is.na(species)) %>% subset_taxa(species != "Homo sapiens")

# Add abundance
k.ps.species@sam_data$Bracken.readct.total <- sample_sums(k.ps.species)

# Plot ordination
ord <- ordinate(transform(k.ps.species, "clr"), method = "MDS", distance = "euclidean")

plot_ordination(k.ps.species, ordination = ord, color = "Sample.Category") +
  geom_point(size = 3, alpha = 1, shape = 21, colour = "black", aes(fill = Sample.Category)) +
  scale_fill_manual(values = sample_palette) +
  theme_classic()

ggsave(file.path("..", "R_output", "k.ps.species_pcoa.png"), width=10, height=8)

# FILTERING PIPELINE ------------------------------------------------------

# FILTERING TAXA:
## DECONTAM
### Not removing taxa with decontam

## ABUNDANCE BASED FILTERING
### Bracken abundance @ 0.01% relative abundance
af.thresh <- c('0.01' = 1e-04)

### this variable is only used for naming right now. the level of filtering is hardcoded like this:
### c("soil", "Ext Blank", "LP Blank") & Copies_in_pool >= max

# FILTERING OUT SAMPLES:
# FEAST
### min oral proportion
min_oral_prop_feast <- 0.30 # 30%
min_oral_prop_decom <- 75

## CUMULATIVE DECAY CURVES
### hard filter
### before AND after decontamination has been performed

# FILTERING BASED ON ABUNDANCE IN ENV CONTROLS
## IN SUMMARY:
### 1) compare relative abundance for each taxa, for each sample, against the relative abundance of matching taxa in controls
### 2) if in any samples, this taxa has a higher relative abundance in controls compared to samples, flag this taxa to be removed from the dataset
### 3) Quality assurance:
#### look through these results to see if this method removes typical oral taxa
#### look to see if typical contaminant species are also removed by this method

# create dummy df to collect results
raw.nums <- tibble('Number of species' = ntaxa(k.ps.species), 'Number of Samples' = nsamples(k.ps.species))
write.csv(raw.nums, file = file.path("..", "R_output", "raw.numbers.csv"), row.names = FALSE)

# RUN FEAST ---------------------------------------------------------------
# important to define what the proportions are here... proportion of taxa present in the original dataset? otherwise 

# prepare otu table
otu_sources_otu <- otu_sources %>% 
  column_to_rownames("TAXID") %>% 
  otu_table(.,taxa_are_rows = T) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) # filter out expty taxa

# add sample data
otu_sources_sam <- bind_cols("SampleID" = sample_names(otu_sources_otu), "SourceSink"="Source") %>%
  #separate(SampleID,into=c("Sample","Env"), sep="\\_", extra="merge", remove=FALSE) %>%
  separate(SampleID,into=c("Sample","Source"), sep="\\_", extra="merge", remove=FALSE) %>%
  mutate(Env = case_when(Source %in% c("human_calculus", "human_plaque", "human_caries") ~ "oral",
                         TRUE ~ Source)) %>%
  column_to_rownames("SampleID") %>% 
  data.frame() %>%  sample_data()

# make a phyloseq object
otu_sources.ps <- phyloseq(otu_sources_otu, otu_sources_sam)
sample_names(otu_sources.ps) <- str_remove(sample_names(otu_sources.ps), "_.*")

# merge them together
SourceSink <- merge_phyloseq(k.ps.species, otu_sources.ps)

# make the accompanying df
SourceSink.df <- SourceSink@sam_data %>% data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  dplyr::select(SampleID, SourceSink, Env) %>% 
  mutate(SourceSink=ifelse(is.na(SourceSink),"Sink", SourceSink)) %>% 
  column_to_rownames("SampleID") %>%  # ROWS MUST BE NAMED AND MATCH WITH SAMPLEID IN OTU TABLE
  arrange(SourceSink) %>%
  mutate(id=ifelse(SourceSink=="Source", NA, row_number()))

SourceSink.otu <- SourceSink %>% abundances() %>% data.frame() %>% mutate_all(as.numeric)

# save tables
feast_dir = file.path("..", "R_output", "feast")
if (!dir.exists(feast_dir)) { dir.create(feast_dir) }

write.table(SourceSink.df, file = file.path(feast_dir, "metadata.feast.txt"),
            row.names = TRUE, col.names = TRUE, na = "NA", quote = FALSE, sep="\t")

write.table(SourceSink.otu, file = file.path(feast_dir, "otu.feast.txt"), quote = FALSE, sep = "\t")

metadata <- Load_metadata(file.path(feast_dir, "metadata.feast.txt"))
otus <- Load_CountMatrix(file.path(feast_dir, "otu.feast.txt"))

if (!file.exists(file.path(feast_dir, "feast_output_source_contributions_matrix.txt"))) {
  FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = feast_dir, outfile="feast_output")
  setwd(wd)
}
FEAST_output <- read.table(file.path(feast_dir, "feast_output_source_contributions_matrix.txt"), quote = "\"", row.names = 1)

# clean up the FEAST output
FEAST_output_tidy <- 
  FEAST_output %>%
  mutate_all(as.numeric) %>%
  rownames_to_column("sample") %>%
  # Simplify sample names
  mutate(sample =  sub("_[^_]+$", "", sample)) %>% 
  # Get source category
  pivot_longer(-sample,names_to = "source",values_to = "proportion") %>%
  # Simplify source names
  mutate(source =  sub("_.*", "", source)) %>% 
  # Add env info and simplify (group calculus and plaque as oral)
  left_join(dplyr::select(rownames_to_column(metadata, "source"), c(source, Env))) %>%
  #mutate(Env = case_when(Env %in% c("human_calculus", "human_plaque") ~ "oral",
  #                       TRUE ~ Env)) %>%
  mutate(Env=replace_na(Env, "unknown") %>% str_remove("human_")) %>%
  # Get 0 values were NAs
  mutate(proportion=as.numeric(proportion),
         proportion=replace_na(proportion,0)) %>%
  # Get sample type
  left_join(dplyr::select(meta, c(LP_ID, Sample.Category, Necropolis)), by=c("sample"="LP_ID")) %>%
  # Categories for plotting
  mutate(FEAST_category = factor(case_when(Sample.Category == "Human" ~ Necropolis,
                                           grepl("Blank", Sample.Category) ~ "Blank",
                                           TRUE ~ Sample.Category),
                                 levels = c("Isola Sacra", "Lucus Feroniae", "Selvicciola", "Blank", "Soil"))) %>%
  # Summarise proportion by source type
  group_by(sample, FEAST_category, Env) %>%
  summarise(proportion = sum(proportion)) %>% ungroup(Env) %>%
  mutate(oral_prop=sum(proportion[Env=="oral"]))  %>%
  # Get source a factor
  mutate(Env = factor(Env, levels = c("oral", "gut", "skin", "soil", "labcontam", "unknown"))) %>%
  ungroup

# Save
write.table(FEAST_output_tidy, file.path(feast_dir, "feast_output_tidy.txt"), quote = FALSE, row.names = FALSE)

# Calculate soil to oral ratio
FEAST_wide <- FEAST_output_tidy %>% dplyr::select(sample, FEAST_category, Env, proportion) %>% 
  group_by(sample, FEAST_category, Env) %>% summarise(proportion = sum(proportion)) %>%
  pivot_wider(names_from = "Env", values_from = "proportion") %>%
  mutate(label = ifelse(oral < min_oral_prop_feast & !FEAST_category %in% c("Soil", "Blank"), "*", NA)) %>%
  mutate(oral_soil_ratio = (oral + 0.01) / (soil + 0.01))  %>% ungroup

write.csv(FEAST_wide,  file.path(feast_dir, "feast_wide.csv"), quote = FALSE, row.names = FALSE)

# rearrange for plotting
FEAST_output_tidy <- FEAST_output_tidy %>%
  mutate(Env = factor(Env, levels = rev(levels(FEAST_output_tidy$Env)))) %>%
  mutate(sample = reorder(sample, oral_prop))

# Plot source composition
p1 <- 
  FEAST_output_tidy %>%
  ggplot(aes(x = sample, y = proportion, fill = Env)) +
  geom_bar(stat = "identity") +
  facet_grid(~FEAST_category, space = "free", scales = "free")+
  # geom_hline(yintercept = 0.95)+
  scale_fill_manual(values = source_palette)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = min_oral_prop_feast) +
  guides(fill = guide_legend(reverse=TRUE))

# Plot
p2 <- 
  FEAST_wide %>% 
  ggplot(aes(x = reorder(sample, oral), y = oral_soil_ratio, fill = oral_soil_ratio)) +
  geom_bar(stat = "identity") +
  # Add label for samples that were removed
  geom_text(aes(label = label, y = 10^-2), size = 3, inherit.aes = TRUE) +
  scale_y_continuous(transform = "log10", labels = label_comma(accuracy = 1)) +
  facet_grid(~FEAST_category, space = "free", scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
                          strip.text = element_blank(),
                          legend.text = element_text(size=10),
                          legend.title = element_blank()) +
  scale_fill_gradientn(colors = colorRampPalette(c("#FF7F0E","#1F77B4"))(100),
                       transform = "log10") +
  xlab("Sample") + ylab("ratio oral/soil")

p_grid <- plot_grid(p1 + labs(tag="A."), p2 + labs(tag="B."), ncol = 1, align = "v", axis = "lr", rel_heights = c(3,2))

ggsave(plot=p_grid, filename = file.path(feast_dir, "feast_composition.png"), width = 20, height = 10, units = "cm")

# PLOT DecOM -----------------------------------

decom_dir = file.path("..", "R_output", "decOM")
if (!dir.exists(decom_dir)) { dir.create(decom_dir) }

decom <- decom %>% select(Sink, starts_with("p_")) %>%
  mutate(across(starts_with("p_"), ~ ifelse(is.na(.x), 0, .x))) %>%
  left_join(dplyr::select(meta, c(LP_ID, Sample.Category, Necropolis)), by=c("Sink"="LP_ID")) %>%
  # Categories for plotting
  mutate(category = factor(case_when(Sample.Category == "Human" ~ Necropolis,
                                     grepl("Blank", Sample.Category) ~ "Blank",
                                     TRUE ~ Sample.Category),
                           levels = c("Isola Sacra", "Lucus Feroniae", "Selvicciola", "Blank", "Soil"))) %>%
  mutate(oral_contam_ratio = (p_aOral + p_mOral + 1)/(p_Skin + p_Sediment.Soil + 1)) %>%
  mutate(Sink = reorder(Sink, (p_mOral + p_aOral)))

decom_long <-  decom %>%
  mutate(label = ifelse((p_aOral + p_mOral) < min_oral_prop_decom & !category %in% c("Soil", "Blank"), "*", NA)) %>%
  pivot_longer(p_Sediment.Soil:p_Unknown, names_to = "Env", values_to = "proportion") %>%
  mutate(Env = str_remove(Env, "p_")) %>%
  # Keep only samples in the OTU table
  filter(Sink %in% sample_names(k.ps.species)) %>%
  mutate(Env = factor(Env, levels = rev(c("aOral", "mOral", "Skin", "Sediment.Soil", "Unknown"))))

write.csv(decom_long,  file.path(decom_dir, "decom_long.csv"), quote = FALSE, row.names = FALSE)

p_decom <- 
  ggplot(data = decom_long, aes(x = Sink, y = proportion, fill = Env)) +
  geom_bar(stat = "identity") +
  # Add label for samples that were removed
  geom_text(aes(label = label, y = 10^-2), size = 3, inherit.aes = TRUE) +
  facet_grid(~category, space = "free", scales = "free") +
  scale_fill_manual(values = source_palette2)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = min_oral_prop_decom) +
  guides(fill = guide_legend(reverse=TRUE))

# Plot
p_decom2 <- 
  decom_long %>% 
  ggplot(aes(x = Sink, y = oral_contam_ratio, fill = oral_contam_ratio)) +
  geom_bar(stat = "identity") +
  # Add label for samples that were removed
  geom_text(aes(label = label, y = 10^-7), size = 3, inherit.aes = TRUE) +
  scale_y_continuous(transform = "log10", labels = label_scientific(accuracy = 1)) +
  facet_grid(~category, space = "free", scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
                          strip.text = element_blank(),
                          legend.text = element_text(size=10),
                          legend.title = element_blank()) +
  scale_fill_gradientn(colors = colorRampPalette(c("#FF7F0E","#1F77B4"))(100),
                       transform = "log10") +
  xlab("Sample") + ylab("ratio oral/contaminant")

p_grid <- plot_grid(p_decom + labs(tag="A."), p_decom2 + labs(tag="B."), ncol = 1, align = "v", axis = "lr", rel_heights = c(3,2))

ggsave(plot=p_grid, filename = file.path(decom_dir, "decom_composition.png"), width = 20, height = 10, units = "cm")

# Get bad samples according to FEAST
assess_samples <-  
  # Get some metadata about samples
  data.frame(k.ps.species@sam_data) %>% dplyr::select(Sample.Category, Necropolis, Tooth.health) %>% rownames_to_column("sample") %>%
  # Add relevant FEAST output
  left_join(dplyr::select(decom, c(Sink, starts_with("p_"), oral_contam_ratio)), by = c("sample" = "Sink")) %>%
  mutate(high_oral_prop = ((p_mOral + p_aOral) > min_oral_prop_decom)) %>%
  mutate(high_oral_prop = factor(case_when(is.na(high_oral_prop) ~ FALSE, TRUE ~ high_oral_prop), levels = c(TRUE, FALSE)))

# Add feast output to phyloseq
new_sam_data <- rownames_to_column(data.frame(k.ps.species@sam_data), "sample") %>%
        left_join(dplyr::select(decom, c(Sink, starts_with("p_"), oral_contam_ratio)), by = c("sample" = "Sink"))

k.ps.species@sam_data <- sample_data(column_to_rownames(new_sam_data, "sample"))
k.ps.species@sam_data$high_oral_prop <- assess_samples$high_oral_prop[match(sample_names(k.ps.species), assess_samples$sample)]

# RUN CUPERDEC ------------------------------------------------------------

taxtable <- k.ps.species@otu_table %>% data.frame() %>% 
  rownames_to_column("Taxon") %>% 
  pivot_longer(cols = -Taxon,names_to = "Sample",values_to = "Count") %>% 
  group_by(Taxon, Sample)

human_high_soil <- assess_samples %>% filter(high_oral_prop == FALSE & Sample.Category == "Human") %>% pull(sample)

metadata_map <- data.frame(k.ps.species@sam_data) %>% rownames_to_column("Sample") %>%
  dplyr::select(Sample,Sample.Category) %>% rename(Sample_Source = Sample.Category) %>%
  # Incorporate feast results: bad samples sould be indicated in Sample_Source column
  mutate(Sample_Source = case_when(Sample %in% human_high_soil ~ "Human_high_soil", 
                                   TRUE ~ Sample_Source))

# Get database for cuperdec
database <- homd.df %>%
  rename(Taxon=TAXID) %>%
  dplyr::select(Taxon) %>%
  mutate(Taxon=as.character(Taxon),
         isolation_source="oral",
         Isolation_Source=ifelse(isolation_source=="oral", TRUE, FALSE)) %>%  
  dplyr::select(-isolation_source)

# save tables
cuperdec_dir = file.path("..", "R_output", "cuperdec")
if (!dir.exists(cuperdec_dir)) { dir.create(cuperdec_dir) }

write.table(taxtable, file = file.path(cuperdec_dir, "taxa_table.txt"),
            row.names = TRUE, col.names = TRUE, na = "NA", quote = FALSE, sep="\t")

write.table(metadata_map, file = file.path(cuperdec_dir, "metadata_map.txt"),
            row.names = TRUE, col.names = TRUE, na = "NA", quote = FALSE, sep="\t")

write.table(database, file = file.path(cuperdec_dir, "database.txt"),
            row.names = TRUE, col.names = TRUE, na = "NA", quote = FALSE, sep="\t")

# Run cuperdec
curves <- calculate_curve(taxtable, database = database) %>%
  print()

filter_result <- simple_filter(curves, percent_threshold = min_oral_prop_decom) %>% print()

plot_cuperdec(curves, metadata_map, filter_result, restrict_x = 200) + theme_classic()
ggsave(file.path(cuperdec_dir, "cuperdec_curves.png"))

# Add which samples didn't pass the cuperdec filter
assess_samples <- assess_samples %>% left_join(filter_result, by = c("sample" = "Sample"))
assess_samples <- rename(assess_samples, c("passed_cuperdec"="Passed"))

# Add info to phyloseq
k.ps.species@sam_data <- cbind(data.frame(k.ps.species@sam_data), rename(column_to_rownames(filter_result, "Sample"), "passed_cuperdec" = "Passed")) %>% sample_data()
k.ps.species@sam_data$passed_cuperdec <- factor(k.ps.species@sam_data$passed_cuperdec)

# Plot concentration to feast and cuperdec output
ggplot(data = data.frame(k.ps.species@sam_data), aes(y = oral_contam_ratio + 0.01, x = Sample.Category, colour = passed_cuperdec)) +
  geom_violin() +
  geom_point(alpha = 0.4, size = 2) +
  scale_y_continuous(trans = "log10") +
  #scale_x_continuous(trans = "log10") +
  theme_classic()

# HEATMAPS ----------------------------

# Colour gradient for abudances
grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4" , "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

# All samples/blanks
dev.off()
png(filename = file.path("..", "R_output", "heatmap_top50_all.png"), width = 40, height = 20, units = "cm", res = 300)
microbiomeutilities::plot_taxa_heatmap(k.ps.species, subset.top=50, transformation="clr",
                                       VariableA=c("Sample.Category", "Necropolis", "high_oral_prop"),
                                       annotation_colors = list("Sample.Category" = sample_palette, "Necropolis" = necr_palette, "high_oral_prop" = c(`TRUE` = "darkgreen", `FALSE` = "white")),
                                       heatcolors = grad_palette)
dev.off()

png(filename = file.path("..", "R_output", "heatmap_top200_all.png"), width = 40, height = 20, units = "cm", res = 300)
microbiomeutilities::plot_taxa_heatmap(k.ps.species, subset.top=200, transformation="clr",
                                       VariableA=c("Sample.Category", "high_oral_prop"),
                                       annotation_colors = list("Sample.Category" = sample_palette, "Necropolis" = necr_palette, "high_oral_prop" = c(`TRUE` = "darkgreen", `FALSE` = "white")),
                                       show_rownames = FALSE,
                                       heatcolors = grad_palette)
dev.off()

# Only human
png(filename = file.path("..", "R_output", "heatmap_top50_human.png"), width = 40, height = 20, units = "cm", res = 300)
microbiomeutilities::plot_taxa_heatmap(subset_samples(k.ps.species, Sample.Category=="Human"), subset.top=50, transformation="clr",
                                       VariableA=c("Necropolis", "Tooth.health"),
                                       annotation_colors = list("Necropolis" = necr_palette, "Tooth.health" = health_palette),
                                       heatcolors = grad_palette)
dev.off()

# RUN DECONTAM ------------------------------------------------------------
## with the following parameters:
### neg : prevalence-based testing. TRUE if sample is a negative control, and FALSE if not (NA entries are not included in the testing). Extraction controls give the best results. If seqtab was provided as a phyloseq object, the name of the appropriate sample-variable in that phyloseq object can be provided.
### method = "either": Contaminants are called if identified by either the frequency or prevalance methods
### conc :  frequency-based testing. A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing. All values must be greater than zero. Zero is assumed to represent the complete absence of DNA. If seqtab was prodivded as a phyloseq object, the name of the appropriate sample-variable in that phyloseq object can be provided. MUST BE GREATER THAN 0!
### threshold : A length-two vector can be provided when using the either or both methods: the first value is the threshold for the frequency test and the second for the prevalence test
### batch	(Optional). factor, or any type coercible to a factor. Default NULL. If provided, should be a vector of length equal to the number of input samples which specifies which batch each sample belongs to (eg. sequencing run). Contaminants identification will be performed independently within each batch. If seqtab was provided as a phyloseq object, the name of the appropriate sample-variable in that phyloseq object can be provided.
### batch.combine	(Optional). Default "minimum". For each input sequence variant (or OTU) the probabilities calculated in each batch are combined into a single probability that is compared to 'codethreshold' to classify contaminants. Valid values: "minimum", "product", "fisher".

# Create folder to save output
decontam_dir = file.path("..", "R_output", "decontamination")
if (!dir.exists(decontam_dir)) { dir.create(decontam_dir) }

# Remove soil and bad samples according to cuperdec
k.ps.species.nosoil <- subset_samples(k.ps.species, Sample.Category != "Soil" & high_oral_prop == TRUE)

negs <- ifelse(k.ps.species.nosoil@sam_data$Sample.Category %in% c("Ext Blank","LP Blank"), TRUE, FALSE)
concs <- k.ps.species.nosoil@sam_data$Library.Preparation.qPCR.Copies.ul. %>% str_remove(" .*") %>%
  as.numeric()

# Write function to test different decontam thresholds
test_decontam <- function(ps, threshold) {
  # Run decontam
  print(paste("Running decontam for threshold", threshold))
  contamdf <- isNotContaminant(ps, neg = negs,
                               method = "prevalence", normalize = TRUE,
                               detailed = TRUE,
                               threshold = threshold)

  # Combine with taxonomic data
  contamdf <- cbind(contamdf, data.frame(ps@tax_table))
  
  # Check if these taxa are present in HOMD and core hominid microbiome
  contamdf$in.homd <- (contamdf$genus %in% homd.df$Genus)
  contamdf$in.core.hominid <- (contamdf$genus %in% str_remove(hominid.df$Taxon, " .*"))
  contamdf$in.contam.list <- (contamdf$genus %in% contaminant_list$genus)
  contamdf <- contamdf %>% mutate(taxon_category =
                                    case_when((in.homd | in.core.hominid) & !in.contam.list ~ "In oral list",
                                              !(in.homd | in.core.hominid) & in.contam.list ~ "In contam list",
                                              (in.homd | in.core.hominid) & in.contam.list ~ "In both",
                                              TRUE ~ "None")) %>%
    mutate(taxon_category = factor(taxon_category, levels = c("In oral list", "In both", "None", "In contam list")))
  contamdf$threshold <- threshold
  # Add taxon name as a different row
  contamdf <- contamdf %>% rownames_to_column("taxon")
  
  # Return table
  return(contamdf)
}

# Iterate over several thresholds
contamdf_all <- data.frame()

# The null hypothesis is that a taxon is a contaminant
# So if a taxon has a p-value below the threshold, it is considered NOT a contaminant
# The higher the threshold, the more taxa pass (the least the contaminants)
for (thr in c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9)) {
  contamdf_all <- rbind(contamdf_all, test_decontam(k.ps.species.nosoil, thr))
}

plot_labs <- contamdf_all %>% group_by(not.contaminant, threshold) %>%
  summarise(count = n_distinct(taxon)) %>%
  mutate(taxon_category = NA)

# Plot
ggplot(data = contamdf_all, aes(x = not.contaminant, fill = taxon_category)) +
  geom_bar(position = "fill") +
  facet_wrap(~threshold) +
  scale_fill_manual(values = c("In oral list" = "#6BBD22", "In both" = "#D4AC27", "None" = "#E2E2E2", "In contam list" = "#C8243E")) +
  theme_classic() +
  geom_label(aes(label = count, y = 0.5), data = plot_labs)

ggsave(file = file.path(decontam_dir, "decontam_tests.png"))

# Write table
write.csv(contamdf_all, file = file.path(decontam_dir, "decontam_results.csv"))

# Not removing any contaminant taxa at this stage
# To assess likely contaminants, get highest threshold for which a taxon is still considered a contaminant
# The higher this value, the stronger the evidence it is a contaminant, so I will just call it p_contaminant although this may not be entirely accurate
assess_taxa <- contamdf_all %>% ungroup %>%
  group_by(taxon, not.contaminant) %>% slice_max(threshold) %>% ungroup() %>% 
  group_by(taxon) %>% slice_min(not.contaminant) %>%
  dplyr::select(taxon, species, starts_with("in."), not.contaminant, threshold)

# RELATIVE ABUNDANCE FILTERING -----------------------------------------------------
# join the orginal phyloseq object with one that has been transformed to relative abundance 
# then censor abundances with relative abundance lower than the threshold on a per-sample basis

# Remove bad samples but keep blanks and soil
k.ps.filtered <- k.ps.species %>% subset_samples(Sample.Category != "Human" | high_oral_prop == TRUE)

# Get relative abundance
compositional_otu <- transform(k.ps.filtered, "compositional") %>% otu_table()

# Set abundances less than the threshold to zero
k.ps.filtered@otu_table[which(compositional_otu < af.thresh)] <- 0

# ABSOLUTE ABUNDANCE FILTERING ------------------------
# get the abundance estimates for controls
bracken.blank <- k.ps.filtered@sam_data %>% data.frame() %>%
  filter(grepl("EB|LP",LP_ID)) %>%
  dplyr::group_by(Extraction.Batch) %>%
  summarise(
  #mean = mean(Bracken.readct.total, na.rm = T),
  max = max(Bracken.readct.total, na.rm = T),
  min = min(Bracken.readct.total, na.rm = T)
  #median = median(Bracken.readct.total, na.rm = T)
  ) %>%
  filter(!is.na(Extraction.Batch))

# filter out samples that were lower in depth compared to max of abundance in controls, by extraction batch
samples.lowseqdepth <- k.ps.filtered@sam_data %>% data.frame() %>%
  left_join(bracken.blank, by = "Extraction.Batch") %>%
  filter(!Sample.Category %in% c("Soil", "Ext Blank", "LP Blank") & Bracken.readct.total <= max) %>%
  pull(LP_ID)

print(samples.lowseqdepth) # No samples have less total abundance than the blanks

# Make sure there are no empty samples or blanks
k.ps.filtered <- subset_taxa(k.ps.filtered, taxa_sums(k.ps.filtered) > 0)
k.ps.filtered <- subset_samples(k.ps.filtered, sample_sums(k.ps.filtered) > 0)

filtered.nums <- tibble('Relative abundance threshold' = names(af.thresh),
                        'Number of Samples Before' = nsamples(k.ps.species),
                        'Number of Taxa Before' = ntaxa(k.ps.species),
                        'Number of Samples After' = nsamples(k.ps.filtered), 
                        'Number of Taxa After' = ntaxa(k.ps.filtered))

write.csv(filtered.nums, file = file.path("..", "R_output", "af.numbers.csv"), row.names = FALSE)

# Specify which taxa were removed during the abundance filtering step
assess_taxa$sufficient_abundance <- assess_taxa$taxon %in% taxa_names(k.ps.filtered)

# Plot ordination
ord <- ordinate(transform(k.ps.filtered, "clr"), method = "MDS", distance = "euclidean")

plot_ordination(k.ps.filtered, ordination = ord, color = "Sample.Category") +
  geom_point(size = 3, alpha = 1, shape = 21, colour = "black", aes(fill = Sample.Category)) +
  scale_fill_manual(values = sample_palette) +
  theme_classic()

ggsave(file.path("..", "R_output", "k.ps.filtered_pcoa.png"), width=10, height=8)

# FILTERING USING SOIL SAMPLES -----------------------------------------------
# for now, we are going to start with the bracken filtered ps object

## Env controls, kraken data

# Get abundances of human dc and soil samples
soil_sample_relabund <- transform(k.ps.filtered, "compositional") %>% psmelt %>%
  dplyr::select(OTU, LP_ID, Abundance, Sample.Category) %>% unique %>%
  filter(Sample.Category %in% c("Soil", "Human")) %>% # we only care about soil + sample comparison
  group_by(OTU) %>%
  filter(sum(Abundance, na.rm = T) > 0) %>%
  ungroup()

# Get soil and humans in different columns
soil_vs_sample <- rename(filter(soil_sample_relabund, Sample.Category == "Human"), c("LP_ID_Human" = "LP_ID", "Abundance_Human" = "Abundance")) %>%
  full_join(rename(filter(soil_sample_relabund, Sample.Category == "Soil"), c("LP_ID_Soil" = "LP_ID", "Abundance_Soil" = "Abundance")), by = "OTU", relationship = "many-to-many") %>%
  mutate(less.abundant.in.sample=(Abundance_Human/Abundance_Soil) <1) %>% dplyr::select(-starts_with("Sample.Category")) %>%
  # always.more.in.env is TRUE if a taxon is more abundant in a soil samples than all calculus samples
  group_by(OTU, LP_ID_Soil) %>%
  dplyr::summarize(always.more.in.env = all(less.abundant.in.sample == "TRUE")) %>%
  # a taxon is identified as environmental if it is most abundant in any one soil sample
  group_by(OTU) %>%
  dplyr::summarise(environmental = any(always.more.in.env == "TRUE"))
  # dplyr::summarise(environmental = any(always.more.in.env == "TRUE"))

print(table(soil_vs_sample$environmental))

# Add to big table
assess_taxa$not.abundant.in.soil <- !(soil_vs_sample$environmental[match(assess_taxa$taxon, soil_vs_sample$OTU)])

assess_taxa <- assess_taxa %>% mutate(taxon_category = case_when(
  (in.homd & in.core.hominid) & !in.contam.list ~ "In oral list",
  (in.homd & in.core.hominid) & in.contam.list ~ "In both",
  !(in.homd & in.core.hominid) & !in.contam.list ~ "None",
  !(in.homd & in.core.hominid) & in.contam.list ~ "In contam list")) %>%
  mutate(taxon_category = factor(taxon_category, levels = c("In oral list", "In both", "None", "In contam list"))) %>%
  # Also get taxon status
  mutate(Status = case_when(!not.contaminant & !sufficient_abundance ~ "Decontam & Abundance",
                            !not.contaminant & !not.abundant.in.soil ~ "Decontam & Soil",
                            !not.contaminant & threshold < 0.5 ~ "Decontam low certainty",
                            !not.contaminant ~ "Decontam high certainty",
                            !sufficient_abundance ~ "Abundance",
                            !not.abundant.in.soil ~ "Soil",
                            TRUE ~ "Passing")) %>%
  mutate(Status = factor(Status,
                         levels = c("Decontam & Soil", "Soil", "Decontam & Abundance", "Abundance", "Decontam high certainty", "Decontam low certainty", "Passing")))

# Write tables for samples and taxa
write.csv(assess_taxa, quote = FALSE, row.names = FALSE, file = file.path(decontam_dir, "assess_taxa.csv"))
write.csv(assess_samples, quote = FALSE, row.names = FALSE, file = file.path(decontam_dir, "assess_samples.csv"))

# Plot
p1 <- ggplot(data = assess_taxa, aes(x = taxon_category, fill = Status)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Decontam & Soil" = "#9E2607", "Soil" = "#C03D1B",
                               "Decontam & Abundance" = "#7E064C", "Abundance" = "#991562",
                               "Decontam high certainty" = "#9E6407", "Decontam low certainty" = "#C0801B",
                               "Passing" = "#6AC022")) +
  theme_classic() + theme(legend.position = "top")

p2 <- ggplot(data = assess_taxa, aes(x = Status, fill = taxon_category)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("In oral list" = "#6BBD22", "In both" = "#D4AC27", "None" = "#E2E2E2", "In contam list" = "#C8243E")) +
  theme_classic() + theme(legend.position = "top")

p <- plot_grid(p1, p2)

ggsave(p, file = file.path(decontam_dir, "taxa_removal_summary.png"), width = 14, height = 8)

# Remove environmental taxa
taxa_to_keep <- soil_vs_sample$OTU[soil_vs_sample$environmental == FALSE]

# Only remove from the samples
k.ps.clean.soil <- k.ps.filtered %>% subset_samples(Sample.Category == "Soil")
k.ps.clean.other <- k.ps.filtered %>% subset_samples(Sample.Category != "Soil")
k.ps.clean.other <-  prune_taxa(taxa_to_keep, k.ps.clean.other)
k.ps.clean <- merge_phyloseq(k.ps.clean.other, k.ps.clean.soil) %>%
  # Also remove RDC_6_H67 which is a duplicate sample
  subset_samples(LP_ID != "RDC_6_H67")
k.ps.clean <- k.ps.clean %>% subset_taxa(taxa_sums(k.ps.clean) > 0)

# Plot ordination
ord <- ordinate(transform(k.ps.clean, "clr"), method = "MDS", distance = "euclidean")

plot_ordination(k.ps.clean, ordination = ord, color = "Sample.Category") +
  geom_point(size = 3, alpha = 1, shape = 21, colour = "black", aes(fill = Sample.Category)) +
  scale_fill_manual(values = sample_palette) +
  theme_classic()

ggsave(file.path("..", "R_output", "k.ps.clean_pcoa.png"), width=10, height=8)

# SUMMARY FILES -----------------------------------------------------
# for comparing to each list of reference taxa...
full.summary <- data.frame(`taxon_category` = "all",
    `flagged` = sum(soil_vs_sample$environmental==TRUE),
    `retained` = sum(soil_vs_sample$environmental==FALSE))

homd.summary <- soil_vs_sample %>%
  mutate(homd_match = case_when(OTU %in% homd.df$TAXID ~ TRUE, TRUE ~ FALSE)) %>%
  filter(homd_match == TRUE) %>%
  summarise(`taxon_category` = "HOMD",
    `flagged` = uniqueN(OTU[environmental == TRUE]),
    `retained` = uniqueN(OTU[environmental == FALSE]))

full.summary <- rbind(full.summary, homd.summary)

hominid.summary <- soil_vs_sample %>%
  mutate(hominid_match = case_when(OTU %in% hominid.df$TAXID ~ TRUE, TRUE ~ FALSE)) %>%
  filter(hominid_match == TRUE) %>%
  summarise(`taxon_category` = "core_hominid",
            `flagged` = uniqueN(OTU[environmental == TRUE]),
            `retained` = uniqueN(OTU[environmental == FALSE]))

full.summary <- rbind(full.summary, hominid.summary)

contaminant.summary <- soil_vs_sample %>% left_join(dplyr::select(psmelt(k.ps.filtered), c(OTU, genus))) %>%
  mutate(contam_match = case_when(genus %in% contaminant_list$genus ~ TRUE, TRUE ~ FALSE)) %>%
  filter(contam_match == TRUE) %>%
  summarise(`taxon_category` = "common_contams",
            `flagged` = uniqueN(OTU[environmental == TRUE]),
            `retained` = uniqueN(OTU[environmental == FALSE]))

full.summary <- rbind(full.summary, contaminant.summary)

full.summary <- full.summary %>% mutate(percentage_flagged = flagged*100/(flagged+retained))

write.csv(full.summary, quote = FALSE, row.names = FALSE, file = file.path(decontam_dir, "decontamination_summary.csv"))

# SUBSET TO ORAL TAXA ------------------------------------------------------
# take taxa found in homd and hominid list, which are not found in contaminants
oral_taxa <- assess_taxa %>% filter((in.homd | in.core.hominid) & !in.contam.list) %>% 
  pull(taxon)

k.ps.oral <- 
  subset_taxa(k.ps.clean, taxa_names(k.ps.clean) %in% oral_taxa)

k.ps.oral <-
  subset_samples(k.ps.oral, sample_sums(k.ps.oral) > 0)

# Plot ordination
ord <- ordinate(transform(k.ps.oral, "clr"), method = "MDS", distance = "euclidean")

plot_ordination(k.ps.oral, ordination = ord, color = "Sample.Category") +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = sample_palette)

ggsave(file.path("..", "R_output", "k.ps.oral_pcoa.png"), width=10, height=8)

# What proportion of the abundance in the clean phyloseq does this subset phyloseq represent?
abundance_retained <-
  k.ps.clean %>% transform("compositional") %>% psmelt() %>%
  dplyr::select(OTU, Sample, Abundance, Sample.Category) %>%
  mutate(retained = OTU %in% taxa_names(k.ps.oral)) %>%
  group_by(Sample, Sample.Category, retained) %>%
  summarise(abundance = sum(Abundance))

ggplot(abundance_retained, aes(x = Sample, y = abundance, fill = retained)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(Sample.Category), scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(`TRUE` = "#6BBD22", `FALSE` = "transparent"))

ggsave(file = file.path(decontam_dir, "oral_abundance_retained.png"))

# SAVE PHYLOSEQ OBJECTS ------------------------------------------------------
ps_list <- c("k.ps", "k.ps.species", "k.ps.filtered", "k.ps.clean", "k.ps.oral")

for (ps_name in ps_list) {
  ps <- get(ps_name)
  saveRDS(ps, file = file.path("..", "R_output", paste0(ps_name, ".RDS")))
}

