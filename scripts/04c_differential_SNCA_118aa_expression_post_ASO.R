# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(readxl)
library(cowplot)

# Arguments---------------------------------------------------------------

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

Transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

plot_transcript_expression <-
  function(data, transcript, samples, genotype, treatment) {
    
    
    # Filter samples to include
    samples_to_include <- samples[samples$Treatment %in% treatment & samples$Mutation %in% genotype, ]$Sample
    
    Expression_per_transcript <-
      data %>%
      dplyr::filter(isoform == transcript) %>% 
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>% 
      dplyr::select(isoform, ORF_length, ORF_seq, all_of(samples_to_include)) %>%
      pivot_longer(!c(isoform, ORF_seq, ORF_length), 
                   names_to = "Sample", 
                   values_to = "count") %>% 
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     Sample, 
                                     Mutation,
                                     Treatment), 
                       by = c("Sample" = "Sample")) %>% 
      aggregate(data = ., count ~ ORF_length + Treatment + Sample + Mutation, FUN = "sum") %>% 
      dplyr::mutate(ORF_length = paste0(ORF_length, "aa"),
                    Sample_ID = gsub("_UT|_M", "", Sample)) 
    
    final_plot <-
      Expression_per_transcript %>% 
      ggplot(aes(x = factor(Treatment, levels = c("UT", "ASO")),
                 y = count
      )) +
      geom_boxplot(fill = "#add8e6", width = 0.5) +
      geom_point(position = position_dodge(width = 0.5)) +
      geom_line(aes(group = Sample_ID),
                color = "grey") + # Add connecting lines
      stat_compare_means(label.x.npc = "center", label = "p.format", paired = F) +
      labs(x = "", 
           y = "Expression per transcript category") +
      scale_y_continuous(labels = function(y) paste0(y, '%')) +
      facet_grid(cols = vars(factor(Mutation, levels = c("Ctrl", "A53T", "SNCAx3"))),
                 rows = vars(ORF_length),
                 scales = "free") +
      theme_bw() +
      theme(plot.title = element_text(size = 14,
                                      face = "bold",
                                      hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "none",
            axis.title = element_text(size = 12, face = "bold"), 
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 10,
                                       face = "bold"),
            axis.text.y = element_text(face = "bold", size = 10))
    
    return(final_plot)
  }

# Main --------------------------------------------------------------------

expression_by_transcript_plot <- plot_transcript_expression(data = Transcripts,
                                                            transcript = "PB.6.191",
                                                            samples = Samples, 
                                                            genotype = c("Ctrl", "A53T", "SNCAx3"), 
                                                            treatment = c("UT", "ASO"))


# Save data ---------------------------------------------------------------
file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = expression_by_transcript_plot,
    filename = paste0("04c_diff_SNCA_118aa_expression_ASO.", ext),
    path = here::here("results", "figures"),
    width = 8,
    height = 6,
    dpi = 600, 
    bg = "white"
  )
}
