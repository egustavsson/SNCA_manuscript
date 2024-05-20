# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(rtracklayer)
library(readxl)

# Arguments---------------------------------------------------------------

args <-
  list(
    path_to_brain_transcripts = here::here("raw_data", "PB_Brain_23042021", "SNCA_RFC1_classification_processed.txt"),
    path_to_iPSC_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

Brain <-
  read.table(args$path_to_brain_transcripts, header = TRUE, sep="\t")

iPSC <-
  read.table(args$path_to_iPSC_transcripts, header = TRUE, sep="\t")

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

plot_ORF_expression <- function(brain_data, cell_data, samples, genotype, treatment, ORFs) {
  
  # Brain regions
  Brain_regions <- c(
    "Spinal_cord",
    "Dorsal_root_ganglion",
    "Medulla_oblongata",
    "Pons",
    "Caudate_nucleus",
    "Thalamus",
    "Hippocampus",
    "Frontal_lobe",
    "Cerebral_cortex",
    "Temporal_lobe",
    "Corpus_callosum",
    "Cerebellum"
  )
  
  # Process Brain data
  Brain_expression_per_ORF <- 
    dplyr::filter(brain_data, ORF_seq %in% ORFs) %>% 
    dplyr::rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)), dplyr::starts_with('NFLR.Clontech_5p..')) %>% 
    dplyr::select(ORF_length, ORF_seq, dplyr::all_of(Brain_regions)) %>% 
    dplyr::mutate(across(dplyr::all_of(Brain_regions), ~(. / sum(.)) * 100)) %>%
    tidyr::pivot_longer(!c(ORF_seq, ORF_length), names_to = "Sample", values_to = "Relative_expression") %>% 
    aggregate(data = ., Relative_expression ~ ORF_length + Sample + ORF_seq, FUN = "sum") %>% 
    dplyr::mutate(ORF_length = paste0(ORF_length, "aa"),
                  source = "Central nervous system") %>% 
    #dplyr::mutate(Sample = gsub("_", " ", Sample)) %>% 
    dplyr::select(ORF_length, Sample, Relative_expression, source)
  
  # iPSC data
  samples_to_include <- samples[samples$Treatment %in% treatment & samples$Mutation %in% genotype, "Sample"]
  
  iPSC_expression_per_ORF <- 
    dplyr::rename_with(cell_data, ~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)), dplyr::starts_with('NFLR.Clontech_5p..')) %>% 
    dplyr::select(ORF_length, ORF_seq, dplyr::all_of(samples_to_include)) %>% 
    dplyr::mutate(across(dplyr::all_of(samples_to_include), ~(. / sum(.)) * 100)) %>%
    tidyr::pivot_longer(!c(ORF_seq, ORF_length), names_to = "Sample", values_to = "Relative_expression") %>% 
    aggregate(data = ., Relative_expression ~ ORF_length + Sample + ORF_seq, FUN = "sum") %>% 
    aggregate(data = ., Relative_expression ~ ORF_length + ORF_seq, FUN = "mean") %>% 
    dplyr::mutate(ORF_length = paste0(ORF_length, "aa"),
                  source = "mDA neurons", 
                  Sample = "Ctrl") %>% 
    dplyr::select(ORF_length, Sample, Relative_expression, source)
  
  # Plot
  Expression_per_ORF_plot <- 
    dplyr::bind_rows(Brain_expression_per_ORF, iPSC_expression_per_ORF) %>%
    dplyr::mutate(Sample = factor(Sample, levels = c("Ctrl", Brain_regions))) %>% 
    ggplot(aes(x = Sample, 
               y = Relative_expression,
               fill = source)) +
    geom_col(width = 0.7, colour = "black", position = position_dodge(width = 0.8)) +
    labs(title = "SNCA open reading frame expression in iPSC-derived mDA neurons and Central nervous system",
         x = "", 
         y = "Relative expression by ORF (%)") +
    scale_y_continuous(labels = function(x) paste0(x, '%')) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_fill_manual(
      values = c("Central nervous system" = "#8D8D8D", "mDA neurons" = "#A9D2DC")
    ) +
    facet_grid(cols = vars(
      factor(source, levels = c("mDA neurons", "Central nervous system"))
             ), 
               rows = vars(ORF_length), 
               scales = "free", 
               space = "free_x") +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 16), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
          axis.text.y = element_text(size = 12))
  
  return(Expression_per_ORF_plot)
}


# Main --------------------------------------------------------------------

# Get all the ORFs of interest from iPSC data
ORF_of_interest <- iPSC %>% 
  dplyr::select(ORF_length, ORF_seq) %>% 
  na.omit() %>% 
  unique() %>% 
  pull(ORF_seq)


expression_by_ORF_plot <- plot_ORF_expression(brain_data = Brain, 
                                              cell_data = iPSC,
                                              samples = Samples, 
                                              genotype = "Ctrl", 
                                              treatment = "UT",
                                              ORFs = ORF_of_interest)


# Save data ---------------------------------------------------------------
file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = expression_by_ORF_plot,
    filename = paste0("06b_ORF_expression.", ext),
    path = here::here("results", "figures"),
    width = 16,
    height = 14,
    dpi = 600, 
    bg = "white"
  )
}
