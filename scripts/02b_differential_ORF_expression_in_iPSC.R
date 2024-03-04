# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(readxl)
library(rtracklayer)
library(ggtranscript)
library(cowplot)

# Arguments---------------------------------------------------------------

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"),
    path_to_gff = here::here("raw_data", "PB_iPSC_04062021", "SNCA_corrected.gtf.cds.gff"),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

Transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

gff <-
  rtracklayer::import(args$path_to_gff)

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

plot_transcript_expression <-
  function(data, gff, samples, genotype, treatment) {
    
    # Filter samples to include
    samples_to_include <- samples[samples$Treatment == treatment & samples$Mutation %in% genotype, ]$Sample
    
    Expression_per_category <-
      data %>%
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>% 
      dplyr::select(isoform, Isoform_class, ORF_length, ORF_seq, all_of(samples_to_include)) %>%
      pivot_longer(!c(isoform, Isoform_class, ORF_seq, ORF_length), 
                   names_to = "Sample", 
                   values_to = "count") %>% 
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     Sample, 
                                     Mutation), 
                       by = c("Sample" = "Sample")) %>% 
      aggregate(data = ., count ~ ORF_length + Sample + Mutation, FUN = "sum") %>% 
      dplyr::mutate(ORF_length = paste0(ORF_length, "aa")) %>% 
      dplyr::mutate(group = ifelse(Mutation %in% c("Ctrl", "A53T"), "Ctrl x A53T", "Ctrl x SNCAx3"))
    
    # need to add the controls to the Ctrl_SNCAx3 comparison
    Expression_per_category <-
      bind_rows(
        Expression_per_category,
        Expression_per_category %>% dplyr::filter(Mutation == "Ctrl") %>% dplyr::mutate(group = "Ctrl x SNCAx3")
      )
    
      # Plot data
    
    plot_list <- list()
    
    for (i in unique(Expression_per_category$group)) {
      
      data <- Expression_per_category %>% dplyr::filter(group == i)
      unique_genotypes <- setdiff(unique(Expression_per_category$Mutation), "Ctrl")
      
      plot_list[[i]] <- 
        data %>%
        ggplot(aes(x = factor(Mutation, levels = c("Ctrl", unique_genotypes)),
                   y = count)) +
        geom_boxplot(fill = "#add8e6", width = 0.5, outlier.shape = NA) +
        geom_point() +
        stat_compare_means(label.x.npc = "center", label = "p.format", paired = F) +
        labs(title = paste("Expression per SNCA ORF in", i, "iPSC mDA neurons"),
             x = "", 
             y = "Relative expression per transcript") +
        scale_y_continuous(labels = function(x) paste0(x, '%')) +
        facet_wrap(
          vars(ORF_length), 
          scales = "free_y",
          nrow = 1) +
        theme_bw() +
        theme(plot.title = element_text(size = 16,
                                        face = "bold",
                                        hjust = 0.5),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              legend.position = "top",
              legend.title = element_blank(),
              strip.text = element_text(size = 12, face = "bold"),
              axis.title = element_text(size = 16), 
              axis.text.x = element_text(angle = 45, 
                                         hjust=1,
                                         size = 12),
              axis.text.y = element_text(face = "bold", size = 12))
     }
    
    return(plot_grid(plot_list[[1]],
                     plot_list[[2]],
                     ncol = 1, 
                     labels = c("A", "B"), 
                     label_size = 16))
  }


# Main --------------------------------------------------------------------

expression_by_ORF_plot <- plot_transcript_expression(data = Transcripts,
                                                     samples = Samples, 
                                                     gff = gff,
                                                     genotype = c("Ctrl", "A53T", "SNCAx3"), 
                                                     treatment = "UT")


# Save data ---------------------------------------------------------------
file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = expression_by_ORF_plot,
    filename = paste0("02b_diff_expression_ctrl_mutant.", ext),
    path = here::here("results", "figures"),
    width = 16,
    height = 10,
    dpi = 600, 
    bg = "white"
  )
}
