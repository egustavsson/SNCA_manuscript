# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(readxl)

# Arguments---------------------------------------------------------------

Gene <- "SNCA"

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", paste0(Gene, "_classification_filtered.txt")),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

Transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

plot_transcripts_per_gene <-
  function(data, gene, gene_name, samples, genotype, treatment, labelling) {
    
    # Filter samples to include
    samples_to_exclude <- samples[samples$Treatment != treatment | samples$Mutation != genotype, ]$Sample
    
    data <- 
      data %>% 
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>% 
      dplyr::select(!samples_to_exclude)
    
    samples <- samples %>% dplyr::filter(!Sample %in% samples_to_exclude)
       
    # fill colour to use
    fill_colour <- c("Coding known (complete match)" = "#045a8d",
                     "Coding known (alternate 3/5 end)" = "#74add1",
                     "Coding novel" = "#4d9221",
                     "NMD novel" = "#d53e4f",
                     "Non-coding known" = "#b2abd2",
                     "Non-coding novel" = "#d8daeb")
    
    
    ## Plot number of transcripts per category ##
    Transcripts_per_category <-
      data %>% 
      dplyr::select(Isoform_class, 
                    associated_gene) %>% 
      dplyr::mutate(associated_gene = gene_name) %>% 
      dplyr::count(associated_gene, Isoform_class)
    
    # if categories are missing populate df
    if(length(c("Coding known (complete match)",
                "Coding known (alternate 3/5 end)",
                "Coding novel",
                "NMD novel",
                "Non-coding known",
                "Non-coding novel")) != length(Transcripts_per_category$Isoform_class)) {
      
      missing <-
        data.frame(associated_gene = gene_name,
                   Isoform_class = setdiff(c("Coding known (complete match)",
                                             "Coding known (alternate 3/5 end)",
                                             "Coding novel",
                                             "NMD novel",
                                             "Non-coding known",
                                             "Non-coding novel"), 
                                           Transcripts_per_category$Isoform_class),
                   n = 0)
      
      Transcripts_per_category <- bind_rows(Transcripts_per_category, missing)
    }
    
    
    Transcripts_per_category_plot <-
      Transcripts_per_category %>% 
      ggplot(aes(x = factor(Isoform_class, 
                            levels = c("Coding known (complete match)",
                                       "Coding known (alternate 3/5 end)",
                                       "Coding novel",
                                       "NMD novel",
                                       "Non-coding known",
                                       "Non-coding novel")), 
                 y = n, 
                 fill = Isoform_class)) +
      geom_col(show.legend = F, colour = "Black") +
      scale_fill_manual(values = fill_colour) +
      scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                                  "Coding known (complete match)" = "Coding known\n(complete match)")) +
      labs(y = "No. unique transcripts", x = "Transcript category") +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(size = 16), 
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 12),
            axis.text.y = element_text(face = "bold", size = 12))
    
    
    ## Plot expression per transcript ##
    Expression_per_transcript <-
      data %>%
      dplyr::select(isoform, Isoform_class, samples$Sample) %>%
      pivot_longer(!c(isoform, Isoform_class), names_to = "Sample", values_to = "NFLR") %>%
      group_by(isoform, Isoform_class) %>%
      summarise(NFLR_mean = mean(NFLR), NFLR_sd = sd(NFLR), .groups = "keep") %>%
      arrange(desc(NFLR_mean)) %>%
      tibble::rowid_to_column(., "isoform_index")
    
    Expression_per_transcript_plot <-
      Expression_per_transcript %>% 
      ggplot(aes(x=isoform_index, fill = Isoform_class)) +
      geom_bar(aes(y = NFLR_mean), color=NA, size=0.3, width=0.8, stat="identity") +
      geom_errorbar(aes(ymin = NFLR_mean-NFLR_sd, ymax = NFLR_mean+NFLR_sd), width = 0.2) +
      scale_fill_manual(values = fill_colour) +
      scale_y_continuous(name = "Transcript expression\nacross samples (%)") +
      scale_x_continuous(name = "Transcripts ranked by expression") +
      guides(fill=guide_legend(title="Transcript category")) +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(size = 16), 
            axis.text.y = element_text(face = "bold", size = 12),
            axis.text.x = element_text(size = 12),
            axis.line = element_line(color = "black", size=0.4))
    
    ## Plot expression per category ##
    Expression_per_category <-
      data %>%
      dplyr::select(Isoform_class, 
                    associated_gene, 
                    samples$Sample) %>%
      pivot_longer(!c(Isoform_class, associated_gene), 
                   names_to = "Sample", 
                   values_to = "count") %>% 
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     Sample, 
                                     Mutation), 
                       by = c("Sample" = "Sample")) %>% 
      aggregate(count ~ Isoform_class + associated_gene + Mutation + Sample,
                data = .,
                FUN = "sum") %>% 
      group_by(Mutation, Isoform_class) %>%
      summarise(mean_count = mean(count), sd_count = sd(count))
    
    categories_to_check <- c("Coding known (complete match)",
                             "Coding known (alternate 3/5 end)",
                             "Coding novel",
                             "NMD novel",
                             "Non-coding known",
                             "Non-coding novel")
    
    # Loop through each factor and check if it is missing in Expression_per_category
    for (variable in categories_to_check) {
      if (!variable %in% Expression_per_category$Isoform_class) {
        # Create a new row with the missing variable, mean_count = 0, sd_count = 0 for all samples
        new_row <- data.frame(Mutation = rep(samples, each = 1),
                              Isoform_class = variable,
                              mean_count = 0,
                              sd_count = 0)
        # Add the new row to Expression_per_category
        Expression_per_category <- rbind(Expression_per_category, new_row)
      }
    }
    
    Expression_per_category_plot <-
      Expression_per_category %>%
      ggplot(aes(factor(Isoform_class, 
                        levels = c("Coding known (complete match)",
                                   "Coding known (alternate 3/5 end)",
                                   "Coding novel",
                                   "NMD novel",
                                   "Non-coding known",
                                   "Non-coding novel")), 
            y = mean_count, 
            fill = Isoform_class)) +
      geom_col(position = "dodge", color = "black") +
      geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                    position = position_dodge(0.9), width = 0.2) +
      labs(x = "Transcript category", 
           y = "Expression per transcript category") +
      scale_fill_manual(values = fill_colour) +
      scale_y_continuous(labels = function(x) paste0(x, '%')) +
      scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                                  "Coding known (complete match)" = "Coding known\n(complete match)")) +
      theme_bw() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            legend.position = "top",
            legend.title = element_blank(),
            axis.title = element_text(size = 16), 
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 12),
            axis.text.y = element_text(face = "bold", size = 12))
    
    plot <- ggpubr::ggarrange(Transcripts_per_category_plot, Expression_per_transcript_plot, Expression_per_category_plot, 
                              nrow = 1, 
                              common.legend = T, align = "h",
                              labels = labelling, 
                              font.label = list(size = 24),
                              widths = c(1, 1.7, 1))
    
    return(annotate_figure(plot, top = text_grob(paste0(gene_name, " transcripts"), 
                                                 color = "black", face = "bold", size = 20))) # return plots as ggarrange
  }

# Main --------------------------------------------------------------------
transcript_plot <-
  plot_transcripts_per_gene(data = Transcripts, 
                            gene = "ENSG00000288563.1", 
                            gene_name = Gene, 
                            samples = Samples, 
                            genotype = "Ctrl", 
                            treatment = "UT",
                            labelling = c("a", "b", "c"))

# Save data ---------------------------------------------------------------
ggsave(plot = transcript_plot, 
       filename = paste0(Gene, "_transcript_plot.png"), 
       path = here::here("results", Gene), 
       width = 16, 
       height = 8, 
       dpi = 600, 
       bg = "white"
)
