#BiocManager::install("microbiome")
#devtools::install_github("microsud/microbiomeutilities")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("MASS")
#BiocManager::install("ANCOMBC")
#devtools::install_github("cozygene/FEAST")
library(tidyverse) # for data cleaning/manipulation
library(phyloseq) # phyloseq objects!
library(microbiome) # for transforming data
library(microbiomeutilities)
library(vegan) # permanova functions
library(cowplot)
library(pairwiseAdonis) # for pairwise permanova
library(ANCOMBC)
library(MASS)
library(ggplot2)
library(ggthemes)
library(FEAST)
library(vegan)
library(microbiome)
library(cuperdec)
library(scales)

# SET THEME AND PALETTES -----------------

sample_palette <- c("Human" = "#3DA417", "Soil" = "#753200", "Ext Blank" = "#C6C3C1", "LP Blank" = "#F3F2F1")

necr_palette <- c("Lucus Feroniae" = "#0D1F61", "Selvicciola" = "#DA8434", "Isola Sacra" = "#2AAE2A")

health_palette <- c("Healthy" = "#5EABAF", "Healthy in caries mouth" = "#DAA334", "Caries" = "#B51C18")

source_palette <- list("oral"="#CC121B", caries = "#7F0006", "gut"="#A97004", "soil"="#3DC02F", "skin"="#EBAD39", "labcontam"="#32569C", "unknown"="#CBCBCB")

wd <- getwd()

# IMPORT DATA -------------------

## Clean and only oral phyloseq object
k.ps.species <- readRDS(file.path("..", "R_output", "k.ps.species.RDS"))
k.ps.clean <- readRDS(file.path("..", "R_output", "k.ps.clean.RDS"))
k.ps.oral <- readRDS(file.path("..", "R_output", "k.ps.oral.RDS"))
assess_taxa <- read.csv(file.path("..", "R_output", "decontamination", "assess_taxa.csv"))

# SUBSET PHYLOSEQ -----------------------
## Work on ps.clean, but remove samples and soil
ps <- k.ps.clean %>% subset_samples(Sample.Category == "Human")

# Also get separate subsets of healthy and caries teeth
ps.healthy <- ps %>% subset_samples(Tooth.health != "Caries")
ps.caries <- ps %>% subset_samples(Tooth.health == "Caries")

ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps.healthy <- prune_taxa(taxa_sums(ps.healthy) > 0, ps.healthy)
ps.caries <- prune_taxa(taxa_sums(ps.caries) > 0, ps.caries)

ps.blanks <- k.ps.clean %>% subset_samples(Sample.Category %in% c("Ext Blank", "LP Blank"))
ps.blanks <- prune_taxa(taxa_sums(ps.blanks) > 0, ps.blanks)

# SAVE PHYLOSEQ OBJECTS ------------------------------------------------------
ps_list <- c("ps", "ps.healthy", "ps.caries", "ps.blanks")

for (ps_name in ps_list) {
  phy <- get(ps_name)
  saveRDS(phy, file = file.path("..", "R_output", paste0(ps_name, ".RDS")))
}

## IMPORT OTHER LISTS
otu_sources <- read.csv(file.path("..", "R_input", "read_bracken_parsed_xfeast.csv"), header = T,comment.char = "",sep=",") %>%
  filter(TAXID %in% taxa_names(ps))

# Lists of oral taxa
# HOMD
homd.df <- read_delim(file.path("..", "R_input", "HOMD_taxon_table2022-08-12_1660320566.txt"), col_names = T, delim = "\t",skip=1) %>%
  dplyr::rename(TAXID=NCBI_taxon_id) %>% 
  distinct(TAXID, .keep_all = TRUE) %>%
  # Keep only oral taxa
  filter(grepl("Oral$", Body_site))

## Write average rel. abundance per sample type
mean_relabunds <-
  k.ps.species %>% transform("compositional") %>%
  psmelt %>% dplyr::select(Sample, OTU, species, Necropolis, Sample.Category, Abundance) %>%
  mutate(sample_type = case_when(Sample.Category == "Human" ~ Necropolis,
                                 grepl("Blank", Sample.Category) ~ "Blank",
                                 TRUE ~ Sample.Category)) %>%
  # Calculate mean across the entire dataset
  group_by(OTU, species) %>% mutate(overall_mean = mean(Abundance)) %>%
  group_by(sample_type, OTU, species, overall_mean) %>% summarise(mean_abund = mean(Abundance)*100) %>%
  pivot_wider(names_from = sample_type, values_from = mean_abund) %>%
  arrange(desc(overall_mean))

# Keep only taxa where the average in at least one group is  more than 0.01%
mean_relabunds <- 
  mean_relabunds %>% filter(if_any(where(is.numeric), ~ . > 0.01))

write.csv(mean_relabunds, file = file.path("..", "R_output", "mean_abundances.csv"), row.names = FALSE)

# Also save OTU table after decontamination
otu_rel_abund <- k.ps.clean %>% transform("compositional") %>%
  psmelt %>% dplyr::select(OTU, Sample, Sample.Category, Necropolis, species, Abundance) %>%
  arrange(Sample.Category, Necropolis, Sample) %>% unique %>%
  pivot_wider(id_cols = c(OTU, species), names_from = Sample, values_from = Abundance)

# Get average abundance per taxon
mean_abunds <- rowMeans(otu_rel_abund[,3:ncol(otu_rel_abund)])

otu_rel_abund <- otu_rel_abund[order(mean_abunds, decreasing = TRUE),]

write.csv(otu_rel_abund, file = file.path("..", "R_output", "otu_clean_relabund.csv"), row.names = FALSE)

# EXPLORE BLANK CONTENT ------------------

# Get abundance rank of each species in each blank
blank_taxa <- ps.blanks %>% psmelt() %>% filter(Abundance > 0) %>%
  dplyr::select(OTU, Sample, Abundance, Extraction.Batch, Library.Preparation.Batch, species, genus, phylum) %>%
  rename(Blank = Sample) %>% group_by(Blank) %>% arrange(desc(Abundance)) %>% mutate(rank_in_blank = dense_rank(desc(Abundance))) %>%
  left_join(dplyr::select(assess_taxa, species, in.homd, in.core.hominid, in.contam.list), by = "species")

blank_oral <- blank_taxa %>% filter(in.homd) %>% filter(!is.na(Extraction.Batch)) 

ranks_in_samples <- blank_oral %>%
  full_join(
    psmelt(ps) %>% group_by(Sample) %>% mutate(rank_in_sample = dense_rank(desc(Abundance))) %>%
      filter(species %in% blank_oral$species) %>% dplyr::select(species, OTU, Sample, rank_in_sample, Extraction.Batch, Library.Preparation.Batch),
    by = c("OTU", "species", "Extraction.Batch"), relationship = "many-to-many") %>%
  # when an abundance is not from a sample within the relevant batches group as other
  mutate(Extraction.Batch = case_when(is.na(Blank) ~ "other", TRUE ~ as.character(Extraction.Batch)))

write.csv(blank_oral_abundances, file.path("..", "R_output", "oral_taxa_in_blanks.csv"))

ggplot(data = ranks_in_samples, aes(x = paste("batch", Extraction.Batch), y = rank_in_sample)) +
  geom_jitter(alpha = 0.8, width = 0.2) +
  geom_point(data = blank_oral, size = 3, colour = "orange", aes(x = paste("batch", Extraction.Batch), y = rank_in_blank)) +
  facet_wrap(~ species, scales = "free") +
  scale_y_reverse() +
  xlab("Extraction batch") + ylab("Abundance rank") +
  theme_classic()

ggsave(filename = file.path("..", "R_output", "oral_taxa_in_blanks.png"), width = 16, height = 15, units = "cm")

# RUN FEAST ---------------------------------------------------------------
# important to define what the proportions are here... proportion of taxa present in the original dataset? otherwise 
#Import kraken-biom of samples (sinks) and sources

# prepare otu table
otu_sources_otu <- otu_sources %>% 
  column_to_rownames("TAXID") %>% 
  otu_table(.,taxa_are_rows = T) %>% 
  filter_taxa(., function(x) sum(x) > 0, TRUE) # filter out empty taxa

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
SourceSink <- merge_phyloseq(ps, otu_sources.ps)
# k.ps.species to have samples, blank and soil objects instead of just samples

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
feast_dir = file.path("..", "R_output", "feast_after_filt")
if (!dir.exists(feast_dir)) { dir.create(feast_dir) }

write.table(SourceSink.df, file = file.path(feast_dir, "metadata.feast.txt"), row.names = T, col.names = T, na = "NA", quote = F, sep="\t")
write.table(SourceSink.otu,file = file.path(feast_dir, "otu.feast.txt"), quote = F, sep = "\t")

metadata <- Load_metadata(file.path(feast_dir, "metadata.feast.txt"))
otus <- Load_CountMatrix(file.path(feast_dir, "otu.feast.txt")) 

if (file.exists(file.path(feast_dir, "feast_output_source_contributions_matrix.txt"))) {
  FEAST_output <- read.table(file.path(feast_dir, "feast_output_source_contributions_matrix.txt"), quote = "\"", row.names = 1)
} else {
  FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = feast_dir, outfile="feast_output")
  setwd(wd)
  FEAST_output <- read.table(file.path(feast_dir, "feast_output_source_contributions_matrix.txt"), quote = "\"", row.names = 1)
}

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
  left_join(dplyr::select(data.frame(ps@sam_data), c(LP_ID, Sample.Category, Necropolis)), by=c("sample"="LP_ID")) %>%
  # Summarise proportion by source type
  group_by(sample, Necropolis, Env) %>%
  summarise(proportion = sum(proportion)) %>% ungroup(Env) %>%
  mutate(oral_prop=sum(proportion[Env %in% c("oral", "caries")]))  %>%
  # Get source a factor
  mutate(Env = factor(Env, levels = c("oral", "caries", "gut", "skin", "soil", "labcontam", "unknown"))) %>%
  ungroup

# Save
write.table(FEAST_output_tidy, file.path(feast_dir, "feast_output_tidy.txt"), quote = FALSE, row.names = FALSE)

# Plot source composition
p1 <- 
  FEAST_output_tidy %>% 
  ggplot(aes(x = reorder(sample, oral_prop), y = proportion, fill = Env)) +
  geom_bar(stat = "identity") +
  facet_grid(~Necropolis, space = "free", scales = "free")+
  # geom_hline(yintercept = 0.95)+
  scale_y_reverse()+
  scale_fill_manual(values = source_palette)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave(plot=p1, filename = file.path(feast_dir, "feast_composition.png"), width = 20, height = 10, units = "cm")

#HEATMAP ---------------------------

# Colour gradient for abudances
grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4" , "#E7C46E","#C24141"))

ps@sam_data$Bracken.readcount.log10 <- log10(ps@sam_data$Bracken.readct.total)

# All samples
dev.off()
png(filename = file.path("..", "R_output", "heatmap_top50_filtered.png"), width = 40, height = 20, units = "cm", res = 300)
microbiomeutilities::plot_taxa_heatmap(ps, subset.top=50, transformation="clr",
                                       VariableA=c("Necropolis", "Tooth.health", "Bracken.readcount.log10"),
                                       annotation_colors = list("Necropolis" = necr_palette, "Tooth.health" = health_palette),
                                       heatcolors = grad_palette(10))
dev.off()


# OVERVIEW OF DATASETS -------------------

# Dominant taxa per necropolis for healthy samples
dom_healthy <- dominant_taxa(ps.healthy, group = "Necropolis", level = "species")$dominant_overview %>%
  arrange(desc(rel.freq))

write.csv(dom_healthy, file = file.path("..", "R_output", "dominant_species_healthy.csv"), row.names = FALSE)

#  Dominant taxa for caries samples
dom_caries <- dominant_taxa(ps.caries, level = "species", group = "Tooth.health")$dominant_overview %>%
  arrange(desc(rel.freq))

id_to_name <- data.frame(ps.caries@tax_table) %>% filter(rownames(.) %in% dom_caries$dominant_taxa) %>% dplyr::select(species)
dom_caries$dominant_taxa <- id_to_name[dom_caries$dominant_taxa, ]

write.csv(dom_caries, file = file.path("..", "R_output", "dominant_species_caries.csv"), row.names = FALSE)

# Abundance vs Prevalence
set.seed(123)
microbiomeutilities::plot_abund_prev(ps)
ggsave(file.path("..", "R_output", "abundance_prevalence.png"))

# ALPHA DIVERSITY ---------------------------

# Rarefaction curves
p1 <- microbiomeutilities::plot_alpha_rcurve(ps.healthy, group = "Necropolis", subsamples = seq(100, 10000, 500)) +
  scale_fill_manual(values = necr_palette, name = "Cemetery") +
  scale_colour_manual(values = necr_palette, name = "Cemetery") +
  geom_vline(xintercept = 4000, linetype = "dashed") +
  annotate("text", x = 4200, y = 200, label = "Threshold for ANOVA test: 4000 reads/sample", angle = 90, size = 4) +
  labs(title = "Rarefaction curves per Necropolis, only for healthy tooth samples")

p2 <- microbiomeutilities::plot_alpha_rcurve(ps, group = "Tooth.health", subsamples = seq(100, 10000, 500)) +
  scale_fill_manual(values = health_palette) +
  scale_colour_manual(values = health_palette)+
  geom_vline(xintercept = 4000, linetype = "dashed") +
  labs(title = "Rarefaction curves per tooth health status, for all samples")

plot_grid(p1, p2, nrow=2, align = "v")

ggsave(file.path("..", "R_output", "alpha_rarefaction.png"), width = 8, height = 10)

# Observed = n° of taxa ; «Shannon» = evenness
alpha <- estimate_richness(ps, measures=c("Observed","Shannon"))

alpha$Necropolis <- ps@sam_data$Necropolis
alpha$Tooth.health <-  ps@sam_data$Tooth.health
alpha$Bracken.readct.total <- ps@sam_data$Bracken.readct.total

alpha_long <- alpha %>%
  pivot_longer(cols = c(Observed, Shannon), names_to = "Measure", values_to = "Value") %>%
  # Note taxa with less than 4000 reads (excluded from tests)
  mutate(label = case_when(Bracken.readct.total < 4000 ~ "*",
                           TRUE ~ ""))

# Plot alpha diversity
jitter = position_jitter(width = 0.1, seed = 123)
ggplot(alpha_long, aes(x = Tooth.health, y = Value)) +
  #geom_boxplot(aes(fill = Tooth.health), outliers = FALSE) +
  geom_point(aes(colour = Tooth.health), size = 4, alpha = 0.5, 
             position = jitter) +
  # Add label for excluded samples from tests
  geom_text(aes(label = label), 
            position = jitter) +
  facet_grid(cols = vars(Necropolis), rows = vars(Measure), scales = "free", space = "free_x") +
  theme_classic() +
  theme(legend.position="none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_colour_manual(values = health_palette) +
  xlab("Tooth health") + ylab("")

ggsave(file.path("..", "R_output", "alphadiv.png"), width = 8, height = 8)

# ANOVA test: does Necropolis affect observed diversity?
# First remove samples with very low read count

alpha_filt <- filter(alpha, Bracken.readct.total > 4000 & Tooth.health != "Caries")

anova_alpha <- aov(data=alpha_filt,
                   Observed ~ Bracken.readct.total + Tooth.health + Necropolis)

## to check if the data are actually normally distributed
hist(MASS::studres(anova_alpha))

shapiro.test(MASS::studres(anova_alpha)) #to check if the distribution is significantly different from a normal distribution

plot(fitted.values(anova_alpha), residuals(anova_alpha))

car::Anova(anova_alpha, type = "II")

# Do the same with Shannon
anova_shannon <- aov(data=alpha_filt,
                     Shannon ~ Bracken.readct.total + Tooth.health + Necropolis)

## to check if the data are actually normally distributed
hist(MASS::studres(anova_shannon))

shapiro.test(MASS::studres(anova_shannon)) #to check if the distribution is significantly different from a normal distribution

plot(fitted.values(anova_shannon), residuals(anova_shannon))

car::Anova(anova_shannon, type = "II")

# Save sample list
write.csv(rownames(alpha_filt), file.path("..", "R_output", "samples_for_alpha_diversity.csv"), row.names = FALSE, quote = FALSE)

# Compare paired taxa
paired_samples <- data.frame(ps@sam_data) %>% group_by(Museum.ID) %>%
  filter(n_distinct(LP_ID) > 1) %>% pull(LP_ID)

# Also, add info to ps
ps@sam_data$paired_samples <- (sample_names(ps) %in% paired_samples)

ps_paired <- prune_samples(paired_samples, ps)

alpha_paired <- estimate_richness(ps_paired, measures=c("Observed","Shannon"))

alpha_paired$Necropolis <- ps_paired@sam_data$Necropolis
alpha_paired$Tooth.health <-  ps_paired@sam_data$Tooth.health
alpha_paired$Bracken.readct.total <- ps_paired@sam_data$Bracken.readct.total
alpha_paired$Individual <- ps_paired@sam_data$Museum.ID


alpha_paired_l <- alpha_paired %>% pivot_longer(cols = c(Observed, Shannon), names_to = "Measure", values_to = "Value")

ggplot(data = alpha_paired_l, aes(x = Tooth.health, y = Value)) +
  geom_point(aes(colour = Tooth.health, size = 4)) +
  facet_grid(rows = vars(Measure), scales = "free", space = "fixed") +
  geom_line(aes(colour = "black", group = Individual)) +
  theme_bw() +
  theme(legend.position="none") +
  scale_colour_manual(values = health_palette) +
  xlab("Tooth health") + ylab("")

ggsave(filename = file.path("..", "R_output", "alpha_healthy_caries_paired.png"))

# LOOK FOR PATHOGENS -----------------------------
pathogens <- homd.df %>%
  dplyr::select(Genus, Species, Disease) %>%
  # Combine genus and species into one
  mutate(Species = paste(Genus, Species, sep = " "),
         Genus = NULL) %>%
  mutate(Taxon_association =  case_when(
    grepl("endocard|heart", Disease) ~ "Putative pathogenicity (heart)",
    grepl("normal|healthy|common|without disease|abundant in caries free|Not usually attributed to disease", Disease) ~ "Putative pathogenicity (other)",
    Species %in% c("Porphyromonas gingivalis", "Tannerella forsythia", "Treponema denticola") ~ "Periodontitis",
    grepl("periodont", Disease) ~ "Periodontitis",
    grepl("caries", Disease) ~ "Caries",
    !is.na(Disease) ~ "Other",
    TRUE ~ "In HOMD")) %>% unique

  ps_pathogens <- ps %>% transform("compositional") %>% psmelt %>% left_join(pathogens, by = c("species" = "Species"), relationship = "many-to-many")

ggplot(data = ps_pathogens, aes(x = Abundance, y = Sample, fill = Taxon_association)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(Necropolis), space = "free", scales = "free")

# ORDINATION AND PERMANOVA --------------------

# Plot ordination
ord <- ordinate(transform(ps, "clr"), method = "MDS", distance = "euclidean")

plot_ordination(ps, ordination = ord, color = "Necropolis", shape = "Tooth.health") +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = necr_palette, name = "Cemetery") +
  scale_shape_manual(values = c("Healthy" = 16, "Healthy in caries mouth" = 6, "Caries" = 4)) +
  theme_classic() +
  stat_ellipse(aes(group = Necropolis), level = 0.95, type = "t") +
  geom_line(aes(group = Museum.ID)) + labs(tag = "A.")

ggsave(file.path("..", "R_output", "final_ps_pcoa.png"), width = 10, height = 6)

plot_ordination(ps, ordination = ord, color = "Tooth.health", shape = "Necropolis") +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = health_palette) +
  theme_classic() +
  stat_ellipse(aes(group = Tooth.health), level = 0.95, type = "t") +
  geom_line(aes(group = Museum.ID))

# PERMANOVA
set.seed(123)

# Compare healthy teeth across Necropolis
perm_ps <- ps.healthy
abund <- t(otu_table(transform(perm_ps, "clr")))
Reads <- perm_ps@sam_data$Bracken.readct.total
Necropolis <- perm_ps@sam_data$Necropolis
Health <- perm_ps@sam_data$Tooth.health

set.seed(123)
perm1 <- adonis2(abund ~ Reads + Health + Necropolis,
                method = "euclidean", na.action = na.omit, nperm = 1000, by = "margin")
perm1

pairwise.adonis(abund, Necropolis, p.adjust.m = "BH", sim.method = "euclidean")

# Compare healthy and unhealthy teeth in different individuals
perm_ps <- ps %>% subset_samples(Tooth.health %in% c("Healthy", "Caries"))
abund <- t(otu_table(transform(perm_ps, "clr")))
Reads <- perm_ps@sam_data$Bracken.readct.total
Necropolis <- perm_ps@sam_data$Necropolis
Health <- perm_ps@sam_data$Tooth.health

perm2 <- adonis2(abund ~ Reads + Health + Necropolis,
                method = "euclidean", na.action = na.omit, nperm = 100, by = "margin")
perm2

# DIFFERENTIAL ABUNDANCE ANALYSIS -------------------------

# group_var
# used to detect structural zeros
# out_cut
# default is 0.05
# if zero counts make up of less than 5% of the proposed Gaussian mixture model,
# it will be detected as outlier zeros and replaced with NA.
# If you believe your dataset is exempt from erroneous data entries, you can also specify out_cut = 0 to disable this detection.
# zero_cut
# used for filtering non-informative taxa
# with 100 samples and zero_cut = 0.90 (default value), taxa with more than 90 zero entries out of 100 taxa will be discarded
# lib_cut
# used for filtering samples
# with lib_cut = 1000, any samples with library size (total observed counts) less than 1000 will be discarded.
# library size varies a lot across different studies, some may have a lot of samples with library size less than 1000.
# In such cases, sticking with the default value will lose a lot of power.
# If you do not want to filter out any sample based on library sizes, you can simply set lib_cut = 0 to disable this function.

# Run ANCOM-BC for comparing LF and IS
ps.ancom <- ps.healthy
ps.ancom <- ps.ancom %>% subset_taxa(taxa_sums(ps.ancom) > 0)
ps.ancom@sam_data$Necropolis <- factor(ps.ancom@sam_data$Necropolis, levels= c("Lucus Feroniae", "Isola Sacra", "Selvicciola"))

if (file.exists(file.path("..", "R_output", "out.rds"))) {
  out <- readRDS(file = file.path("..", "R_output", "out.rds"))
} else {
  out <- ancombc2(data = ps.ancom,
                 fix_formula = "Necropolis", # Here you add in the form of a formula your variable of interest (Necropolis) and maybe other variables you want to account for
                 tax_level = "species", # You can change this to other levels in your phyloseq. Make sure to spell it how it is in your phyloseq (e.g. species not Species). If you go above family it recommended to change neg_l (below) to TRUE
                 p_adj_method = "BH", prv_cut = 0.10, # Some parameters. prv_cut=0.10 means that taxa found in less than 10% of samples are ignored
                 group="Necropolis",
                 struc_zero = TRUE,
                 lib_cut = 1000,
                 verbose = TRUE)
  
  saveRDS(out, file = file.path("..", "R_output", "out.rds"))
}

## Get results
res <- out$res

write.csv(res, file.path("..", "R_output", "ancom_res.csv"), quote = FALSE, row.names = FALSE)

colnames(res) <- colnames(res) %>% str_remove_all("Necropolis") %>%
  make.names()

# Note which diff. abundant taxa also pass sensitivity tests
res_filt <- res %>% filter(diff_Isola.Sacra == TRUE | diff_Selvicciola == TRUE)

## Get structural zeros table
zero_ind <- out$zero_ind %>%
  filter(taxon %in% res_filt$taxon)

# Give columns somewhat easier names
colnames(zero_ind) <- str_remove(colnames(zero_ind), "structural_zero .* = ") %>%
  str_remove(pattern = "[)]") %>% make.names() 

# Get additional metadata from ps object
ps.res <- ps.ancom %>% transform("compositional") %>%
  subset_taxa(species %in% res_filt$taxon) %>% psmelt %>%
  # Join with 
  left_join(res_filt, by = c("species" = "taxon")) %>%
  dplyr::select(Sample, species, Necropolis, Abundance, contains("diffabund_")) %>% 
  unique

# Plot

box <-
  ggplot(data = ps.res, aes(x = Abundance * 100 + 0.01, y = Necropolis, fill = Necropolis, colour = Necropolis)) +
  geom_boxplot(outliers = FALSE, fill = "white")+
  geom_point(size = 4, alpha = 0.5, position = position_jitter(height = 0.2, width = 0.01)) +
  facet_grid(rows = vars(str_replace(species, " ", "\n")), space = "free", switch = "both") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = "transparent"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12)) +
  scale_x_continuous(transform = "log10") + ylab("") + xlab("% Abundance") +
 # scale_fill_manual(values = sapply(necr_palette, colorspace::lighten, amount = 0.5))+
  scale_colour_manual(values = necr_palette)

# LFC barplot
# Get lfc and standard errors for the two comparisons
lfc <- res_filt %>% dplyr::select(taxon, contains("lfc_"), contains("se_")) %>%
  pivot_longer(contains("lfc_") | contains("se_"), names_to = "Necropolis", values_to = "value") %>%
  separate(Necropolis, into = c("measure", "Necropolis"), sep = "_") %>%
  # Remove intercept lfc
  filter(Necropolis != ".Intercept.") %>%
  pivot_wider(names_from = "measure", values_from = "value") %>%
  mutate(Necropolis = gsub(".", " ", Necropolis, fixed=TRUE))

# Indicate which comparisons passed the sensitivity test
ss <- res_filt %>% 
  # Which taxa are diff abund and pass ss
  mutate(diffpass_Isola.Sacra = diff_Isola.Sacra & passed_ss_Isola.Sacra,
         diffpass_Selvicciola = diff_Selvicciola & passed_ss_Selvicciola) %>%
  dplyr::select(taxon, contains("diffpass_")) %>% 
  pivot_longer(contains("diffpass_"), names_to = "Necropolis", values_to = "diffpass") %>%
  mutate(Necropolis = gsub(".", " ", str_remove(Necropolis, "diffpass_"), fixed = TRUE)) %>%
  filter(!grepl("Intercept", Necropolis)) %>%
  mutate(label = case_when(diffpass ~ "*", TRUE ~ ""))

lfc <- left_join(lfc, ss)

bar <-
  ggplot(data = lfc, aes(x = lfc, y = taxon, fill = Necropolis)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbarh(aes(xmin = lfc - se, 
                      xmax = lfc + se,
                      y = taxon, height = 0.4),
                   position = position_dodge(1)) + ylab("") +
  geom_text(aes(label = label), position = position_dodge(1), size = 8) +
  scale_fill_manual(values = necr_palette) +
  scale_y_discrete(limits = rev(levels(factor(res_filt$taxon)))) +
  theme_classic() + theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          legend.position = "bottom") + 
  geom_vline(xintercept = 0, size = 2) + 
  xlab("Log fold change compared to Lucus Feroniae")

plot_grid(box + labs(tag = "B."),
          bar + labs(tag = "C."),
          align = "h", axis = "tb", ncol = 2, nrow = 1, rel_widths = c(3, 2))

ggsave(file=file.path("..", "R_output", "ancom_plot.png"), width=10, height=7)

ps.res %>% group_by(species) %>% summarise(mean = mean(Abundance[Abundance!=0])*100,
                                                       prevalence=sum(Abundance != 0)/n_distinct(Sample))

# IDENTIFY TAXA TO RUN MAP DAMAGE ON ------------------------------------------

# Get the 20 most abundant taxa per Necropolis
mapdamage_taxa <- transform(ps, "clr") %>% psmelt %>%
  # First get mean abundance per Necropolis
  group_by(Necropolis, species) %>% summarise(Necropolis_mean = mean(Abundance)) %>%
  slice_max(Necropolis_mean, n = 20) %>% pull(species) %>% unique

mapdamage_taxa <- c(mapdamage_taxa, blank_oral$species, res_filt$taxon) %>% unique

# Identify which sample to use for mapping
mapdamage_taxa_samples <- psmelt(ps) %>% filter(species %in% mapdamage_taxa) %>%
  group_by(OTU, species) %>% filter(Abundance == max(Abundance)) %>% dplyr::select(OTU, species, Sample, Abundance) # This excludes Staphylococcus hominis but it was in very low abundances in samples

write.csv(mapdamage_taxa_samples, file.path("..", "R_output", "mapdamage_taxa_samples.csv"), row.names = FALSE, quote = FALSE)

# Save list of samples used for this analysis

write.csv(sample_names(ps), file.path("..", "R_output", "samples_for_taxonomic_analysis.csv"), row.names = FALSE, quote = FALSE)
