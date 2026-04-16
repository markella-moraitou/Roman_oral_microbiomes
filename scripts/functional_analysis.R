#load packages
library(readr)
library(stringr)
library(ape)
library(vegan)
library(phyloseq)
library(ggpubr)
library(tidyverse)
library(microbiome)
library(reshape2)
library(pairwiseAdonis)
library(ANCOMBC)
library(cowplot)

# SET THEME AND PALETTES ------------------------------------

sample_palette <- c("Human" = "#3DA417", "Soil" = "#753200", "Ext Blank" = "#C6C3C1", "LP Blank" = "#F3F2F1")

necr_palette <- c("Lucus Feroniae" = "#0D1F61", "Selvicciola" = "#DA8434", "Isola Sacra" = "#2AAE2A")

health_palette <- c("Healthy" = "#5EABAF", "Healthy in caries mouth" = "#DAA334", "Caries" = "#B51C18")

source_palette <- list("oral"="#CC121B", caries = "#7F0006", "gut"="#A97004", "soil"="#3DC02F", "skin"="#EBAD39", "labcontam"="#32569C", "unknown"="#CBCBCB")

wd <- getwd()

subdir <- normalizePath(file.path("..", "R_output", "functional"))
if (!dir.exists(subdir)) { dir.create(subdir) }

# IMPORT DATA -------------------------------------------------

#Import unstratified GO table
GO_terms <- read_tsv(file.path("..", "R_input", "all_GOterms_cpm_renamed.tsv")) %>% as.data.frame

#Import taxon phyloseq
k.ps.clean <- readRDS(file.path("..", "R_output", "k.ps.clean.RDS"))
ps <- readRDS(file.path("..", "R_output", "ps.RDS"))

#Import taxonomy table
micr_taxonomy <- read.csv(file.path("..", "R_output", "microbial_taxonomy.csv"))

# Import decontamination results
#assess_taxa <- read.csv(file.path("..", "R_output", "decontamination", "assess_taxa.csv"))

# Keep only unstratified GO terms -------------------------------------------------

#Rename columns
colnames(GO_terms) <- colnames(GO_terms) %>% str_remove("_Abundance-RPKs")

# Remove UNMAPPED, UNGROUPED and stratified output
GO_unstr <- 
  GO_terms %>% separate(`# Gene Family`, into = c("GO", "Taxon"), sep = "\\|") %>%
  # Remove unmapped and unstratified data
  filter(GO != "UNMAPPED" & GO != "UNGROUPED" & is.na(Taxon)) %>%
  mutate(Taxon = NULL) %>%
  column_to_rownames("GO")

#Keep only samples that were retained as 'good' samples under community-level taxonomic analysis
GO_unstr <- GO_unstr[,which(colnames(GO_unstr) %in% sample_names(k.ps.clean))]

#Remove empty rows
GO_unstr <- GO_unstr[which(rowSums(GO_unstr)>0),]

### Create phyloseq object ###
GO_phyloseq <- phyloseq(otu_table(GO_unstr, taxa_are_rows = TRUE), sample_data(k.ps.clean))

# Remove empty samples and blanks
GO_phyloseq <- subset_samples(GO_phyloseq, sample_sums(GO_phyloseq) > 0) %>% subset_samples(Sample.Category %in% c("Human", "Soil"))
GO_phyloseq <- prune_taxa(taxa_sums(GO_phyloseq) > 0, GO_phyloseq)

#Create a subset of the phyloseq object that contains only biological processes
GO_BP_phyloseq <- prune_taxa(grepl(" \\[BP\\] ", taxa_names(GO_phyloseq)), GO_phyloseq)

# Remove empty samples
GO_BP_phyloseq <- subset_samples(GO_BP_phyloseq, sample_sums(GO_BP_phyloseq) > 0)
GO_BP_phyloseq <- prune_taxa(taxa_sums(GO_BP_phyloseq) > 0, GO_BP_phyloseq)

#save object
saveRDS(GO_phyloseq, file.path(subdir, "GO_phyloseq"))
saveRDS(GO_BP_phyloseq, file.path(subdir, "GO_BP_phyloseq"))

## CLR normalise

GO_phyloseq_clr <- transform(GO_phyloseq, "clr")
GO_BP_phyloseq_clr <- transform(GO_BP_phyloseq, "clr")

# PCoA ----------------

## Soil and human ##

## All GOs
#Calculate euclidean distance matrix
GO_eucl <- ordinate(GO_phyloseq_clr, method="PCoA", distance="euclidean")

#Plot and save ordinations
plot_ordination(GO_phyloseq_clr, GO_eucl, color = "Sample.Category", shape = "Tooth.health") +
  geom_point(size = 5, alpha = 1) +
  scale_colour_manual(values = sample_palette) +
  scale_shape_manual(values = c("Healthy" = 16, "Healthy in caries mouth" = 6, "Caries" = 4), na.value = 16) +
  stat_ellipse(aes(group = Sample.Category, colour = Sample.Category), level = 0.95, type = "t") +
  geom_vline(xintercept = -30, linetype = "dotted") +
  theme_classic()

ggsave(file.path(subdir, "GO_ordination_all.png"), width=10, height=8)

## Only BP
#Calculate euclidean distance matrix
GO_eucl <- ordinate(GO_BP_phyloseq_clr, method="PCoA", distance="euclidean")

#Plot and save ordinations
plot_ordination(GO_BP_phyloseq_clr, GO_eucl, color = "Sample.Category") +
  geom_point(size = 5, alpha = 1) +
  scale_colour_manual(values = sample_palette) +
  stat_ellipse(aes(group = Sample.Category, colour = Sample.Category), level = 0.95, type = "t") +
  geom_vline(xintercept = -10, linetype = "dotted") +
  theme_classic()

ggsave(file.path(subdir, "BP_ordination_all.png"), width=10, height=8)

## Remove soil and bad samples
# Identify samples that are within the soil ellipse and outside the human ellipse (Axis 1 < -70)
remove_samples <- GO_eucl$vectors %>% as.data.frame %>% filter(Axis.1 < -10) %>% rownames
remove_samples <- data.frame(sample = remove_samples) %>% filter(grepl("_H", sample))

remove_samples$Necropolis <- GO_BP_phyloseq@sam_data$Necropolis[match(remove_samples$sample, sample_names(GO_BP_phyloseq))]

write.csv(remove_samples, file = file.path(subdir, "removed_samples.csv"))

# Remove samples and soil
GO_filt <- GO_phyloseq_clr %>% subset_samples(!sample_names(GO_phyloseq_clr) %in% remove_samples$sample)
GO_filt <- subset_samples(GO_filt, Sample.Category == "Human")

GO_BP_filt <- GO_BP_phyloseq_clr %>% subset_samples(!sample_names(GO_BP_phyloseq_clr) %in% remove_samples$sample)
GO_BP_filt <- subset_samples(GO_BP_filt, Sample.Category == "Human")

## Run ordination with filtered data
#Calculate euclidean distance matrix
GO_eucl <- ordinate(GO_filt, method="PCoA", distance="euclidean")

plot_ordination(GO_filt, GO_eucl, color = "Necropolis", shape = "Tooth.health") +
  geom_point(size = 5) +
  scale_color_manual(values = necr_palette, name = "Cemetery") +
  scale_shape_manual(values = c("Healthy" = 16, "Healthy in caries mouth" = 6, "Caries" = 4)) +
  #scale_colour_viridis_b(trans = "log10") +
  theme_classic() +
  stat_ellipse(aes(group = Necropolis), level = 0.95, type = "t") +
  geom_line(aes(group = Museum.ID))

ggsave(file.path(subdir, "GO_ordination_filt.png"), width=10, height=8)

#Calculate euclidean distance matrix
GO_eucl <- ordinate(GO_BP_filt, method="PCoA", distance="euclidean")

plot_ordination(GO_BP_filt, GO_eucl, color = "Necropolis", shape = "Tooth.health") +
  geom_point(size = 5) +
  scale_color_manual(values = necr_palette, name = "Cemetery") +
  scale_shape_manual(values = c("Healthy" = 16, "Healthy in caries mouth" = 6, "Caries" = 4)) +
  theme_classic() +
  stat_ellipse(aes(group = Necropolis), level = 0.95, type = "t") +
  geom_line(aes(group = Museum.ID)) + labs(tag = "A.")

ggsave(file.path(subdir, "BP_ordination_filt.png"), width=10, height=8)

# Get unstratified biological processes subset
GO_unstr_BP <- GO_unstr[which(grepl("\\[BP\\]", rownames(GO_unstr))), sample_names(GO_BP_filt)] %>%
  rownames_to_column() %>%
  separate(rowname, sep = ": .BP. ", into = c("GO Term ID", "Function"))

write.csv(GO_unstr_BP, file = file.path(subdir, "GO_unstr_BP.csv"), row.names = FALSE, quote = TRUE)

# PERMANOVA ----------------------

# Compare healthy teeth across Necropolis
perm_ps <- GO_BP_filt %>% subset_samples(Tooth.health != "Caries")
abund <- t(otu_table(perm_ps))
Reads <- perm_ps@sam_data$Bracken.readct.total
Necropolis <- perm_ps@sam_data$Necropolis
Health <- perm_ps@sam_data$Tooth.health

set.seed(123)

GO_BP_model <- adonis2(abund ~ Reads + Health + Necropolis,
                 method = "euclidean", na.action = na.omit, nperm = 1000, by = "margin")
GO_BP_model

GO_BP_pair <- pairwise.adonis(abund, Necropolis, p.adjust.m = "BH", sim.method = "euclidean")

GO_BP_pair

# DIFFERENTIAL ABUNDANCE ANALYSIS -------------------------

# Run ANCOM-BC for comparing LF and IS
ps.ancom <- GO_BP_phyloseq %>% subset_samples(sample_names(GO_BP_phyloseq) %in% sample_names(perm_ps))
ps.ancom <- ps.ancom %>% prune_taxa(taxa_sums(ps.ancom) > 0, .)
ps.ancom@sam_data$Necropolis <- factor(ps.ancom@sam_data$Necropolis, levels= c("Lucus Feroniae", "Isola Sacra", "Selvicciola"))

GO_tax_table <- data.frame(rownames = taxa_names(ps.ancom)) %>%
  separate(rownames, into = c("GO_code", "GO_name"), sep = ": \\[BP\\] ", remove = FALSE) %>%
  column_to_rownames("rownames") %>%
  as.matrix %>% tax_table

ps.ancom@tax_table <- GO_tax_table

if (file.exists(file.path(subdir, "GO_out.rds"))) {
  out <- readRDS(file = file.path(subdir, "GO_out.rds"))
} else {
  out <- ancombc2(data = ps.ancom,
                  fix_formula = "Necropolis",
                  tax_level = "GO_name", 
                  p_adj_method = "BH", prv_cut = 0.2, 
                  group="Necropolis",
                  struc_zero = FALSE,
                  lib_cut = 1000,
                  verbose = TRUE)
  
  saveRDS(out, file = file.path(subdir, "GO_out.rds"))
}

## Get results
res <- out$res

colnames(res) <- colnames(res) %>% str_remove_all("Necropolis") %>%
  make.names()

write.csv(res, file.path(subdir, "GO_ancom_res.csv"), quote = FALSE, row.names = FALSE)

# Note which diff. abundant taxa also pass sensitivity tests
res_filt <- res %>% filter((diff_Isola.Sacra & passed_ss_Isola.Sacra) | (diff_Selvicciola & passed_ss_Selvicciola) )

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
  theme_classic() + theme(axis.ticks.y = element_blank(),
                          legend.position = "bottom") + 
  geom_vline(xintercept = 0, size = 2) + 
  xlab("Log fold change compared to Lucus Feroniae")

ggsave(file=file.path(subdir, "ancom_plot.png"), width=8, height=4)
