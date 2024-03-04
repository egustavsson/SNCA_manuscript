# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(rtracklayer)

# Arguments---------------------------------------------------------------

args <-
  list(
    path_to_brain_transcripts = here::here("raw_data", "PB_Brain_23042021", "SNCA_classification_processed.txt"),
    path_to_iPSC_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt")
  )

# Load data ---------------------------------------------------------------

Brain <-
  read.table(args$path_to_brain_transcripts, header = TRUE, sep="\t")

iPSC <-
  read.table(args$path_to_iPSC_transcripts, header = TRUE, sep="\t")

# Functions ---------------------------------------------------------------

plot_ORF_expression <-
  function(data, ORFs) {
    
    filtered_ORFs <- data %>% dplyr::filter(ORF_seq %in% ORFs)
    
    
    Expression_per_ORF <-
      filtered_ORFs %>%
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>% 
      dplyr::select(ORF_length, 
                    ORF_seq, 
                    Caudate_nucleus,
                    Cerebellum,
                    Cerebral_cortex,
                    Corpus_callosum,
                    Dorsal_root_ganglion,
                    Frontal_lobe,
                    Medulla_oblongata,
                    Pons,
                    Spinal_cord,
                    Temporal_lobe,
                    Thalamus) %>%
      pivot_longer(!c(ORF_seq, ORF_length), 
                   names_to = "Sample", 
                   values_to = "Relative_expression") %>% 
      aggregate(data = ., Relative_expression ~ ORF_length + Sample + ORF_seq, FUN = "sum") %>% 
      dplyr::mutate(ORF_length = paste0(ORF_length, "aa"),
                    Sample = gsub("_", " ", Sample))
    
    # Plot
    Expression_per_ORF_plot <- 
      Expression_per_ORF %>%
      ggplot(aes(x = Sample,
                 y = Relative_expression)) +
      geom_col(width = 0.7, colour = "black", position = position_dodge(width = 0.8), fill = "#add8e6") +
      labs(title = "SNCA open reading frame expression in brain",
           x = "", 
           y = "Relative expression by ORF (%)") +
      scale_y_continuous(labels = function(x) paste0(x, '%')) +
      facet_grid(vars(ORF_length), scales = "free_y", space = "free_x") +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "top",
            legend.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16), 
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 12, 
                                       face = "bold"),
            axis.text.y = element_text(size = 12))
    
    
    return(Expression_per_ORF_plot)
  }


# Main --------------------------------------------------------------------

# Get all the ORFs of interest from iPSC data
ORF_of_interest <- iPSC %>% 
  dplyr::select(ORF_length, ORF_seq) %>% 
  na.omit() %>% 
  unique() %>% 
  .$ORF_seq


expression_by_ORF_plot <- plot_ORF_expression(data = Brain,
                                              ORFs = ORF_of_interest)


# Save data ---------------------------------------------------------------
file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = expression_by_ORF_plot,
    filename = paste0("06b_ORF_expression_in_brain.", ext),
    path = here::here("results", "figures"),
    width = 8,
    height = 12,
    dpi = 600, 
    bg = "white"
  )
}
