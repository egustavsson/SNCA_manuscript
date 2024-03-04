# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(readxl)

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

plot_transcript_category_expression <-
  function(data, samples, genotype, treatment) {
    
    fill_colour <- c("Coding known (complete match)" = "#045a8d",
                     "Coding known (alternate 3/5 end)" = "#74add1",
                     "Coding novel" = "#4d9221",
                     "NMD novel" = "#d53e4f",
                     "Non-coding known" = "#b2abd2",
                     "Non-coding novel" = "#d8daeb")
    
    # Filter samples to include
    samples_to_include <- samples[samples$Treatment == treatment & samples$Mutation %in% genotype, ]$Sample
    
    Expression_per_category <-
      data %>%
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>% 
      dplyr::select(Isoform_class, all_of(samples_to_include)) %>%
      pivot_longer(!c(Isoform_class), 
                   names_to = "Sample", 
                   values_to = "count") %>% 
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     Sample, 
                                     Mutation), 
                       by = c("Sample" = "Sample")) %>% 
      aggregate(count ~ Isoform_class + Mutation + Sample,
                data = .,
                FUN = "sum") %>% 
      dplyr::mutate(case_ctrl = ifelse(Mutation == "Ctrl", "Control", "Mutant")) %>% 
      group_by(case_ctrl, Isoform_class)
    
    # Plot data
    Expression_per_category_plot <- 
      Expression_per_category %>%
      ggplot(aes(factor(Mutation,
                        levels = c("Ctrl",
                                   "SNCAx3",
                                   "A53T")),
                 y = count)) +
      geom_boxplot(aes(fill = Isoform_class), width = 0.5, outlier.shape = NA) +
      geom_point() +
      stat_compare_means(label.x.npc = "center") +
      labs(x = "Transcript category", 
           y = "Expression per transcript category") +
      scale_fill_manual(values = fill_colour) +
      scale_y_continuous(labels = function(x) paste0(x, '%')) +
      scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                                  "Coding known (complete match)" = "Coding known\n(complete match)")) +
      facet_wrap(vars(Isoform_class)) +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "top",
            legend.title = element_blank(),
            strip.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16), 
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 12),
            axis.text.y = element_text(face = "bold", size = 12))
    
    
    return(Expression_per_category_plot)
  }
    
# Main --------------------------------------------------------------------

expression_by_category_plot <- plot_transcript_category_expression(data = Transcripts,
                                                                   samples = Samples, 
                                                                   genotype = c("Ctrl", "A53T", "SNCAx3"), 
                                                                   treatment = "UT")


# Save data ---------------------------------------------------------------
file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = expression_by_category_plot,
    filename = paste0("02a_diff_expression_category_ctrl_UT.", ext),
    path = here::here("results", "figures"),
    width = 10,
    height = 8,
    dpi = 600, 
    bg = "white"
  )
}
