# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggforce)
library(readxl)

# Arguments---------------------------------------------------------------

Gene <- "SNCA"

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", paste0(Gene, "_classification_processed.txt")),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

sqanti_class <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

sample_info <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

transcripts_depreciation_curve <-
  
  function(sqanti_output, gene_name, samples) {
    
    # Reduce the data frame to only include samples of interest
    # and calculate NFLR_mean with only the samples included
    sqanti_output <- 
      sqanti_output %>% 
      select(isoform, matches(paste0("^NFLR.*(", paste(samples, collapse = "|"), ")"))) %>% 
      dplyr::mutate(NFLR_mean = rowMeans(across(starts_with("NFLR"))))
    
    counts.long <- sqanti_output %>%
      dplyr::select(isoform, 
                    starts_with("NFLR")) %>%
      pivot_longer(cols = c(-isoform), 
                   names_to = "sample", 
                   values_to = "read_count")
    
    counts.long$sample <- counts.long$sample %>% 
      str_replace("NFLR.*\\.", "") %>% 
      str_replace("_3p", "")
    
    threshs <- seq(from = 0, to = max(sqanti_output$NFLR_mean), by = 0.05)
    
    res <- data.frame()
    
    for(i in threshs){
      transcripts_per_sample <- counts.long %>%
        dplyr::filter(read_count >=i) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(n_transcripts = dplyr::n(), 
                         .groups = "keep") %>%
        dplyr::mutate(NFLR_threshold = i)
      
      res <- rbind(res, transcripts_per_sample)
      
      rm(transcripts_per_sample)
      
    }
    
    df <- res %>%
      dplyr::filter(sample %in% "NFLR_mean")
    
    sd <- res %>%
      dplyr::filter(!sample %in% "NFLR_mean") %>%
      dplyr::group_by(NFLR_threshold) %>%
      dplyr::summarise(sd = sd(n_transcripts), 
                       .groups = "keep")
    
    df <- left_join(df, 
                    sd, 
                    by = "NFLR_threshold")
    
    # Plot data
    plot <-
      df %>% 
      ggplot(aes(x=NFLR_threshold, 
                 y = n_transcripts)) +
      geom_line(group=1) +
      geom_ribbon(aes(NFLR_threshold, 
                      ymax = n_transcripts + sd, 
                      ymin = n_transcripts - sd), 
                  alpha = 0.2, 
                  fill = "blue", group=1) +
      scale_x_continuous(name = expression(paste("Normalised expression (", NFLR[T], ")"))) +
      scale_y_continuous(name = "No. of transcripts detected") +
      guides(fill="none") +
      facet_zoom(xlim = c(0, 5),
                 ylim = c(0, 50),
                 horizontal = F) +
      
      theme(plot.title = element_text(face = "bold",
                                      size = 16,
                                      hjust = 0.5),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size=0.4),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(face = "bold",
                                       size = 10),
            axis.text.y = element_text(face = "bold",
                                       size = 10),
            axis.title.y = element_text(size = 14),
            strip.text.x = element_text(face = "bold",
                                        size = 12),
            legend.title = element_blank())
    
    return(plot)
  }

# Main --------------------------------------------------------------------

# Samples to run plot
samples_to_plot <- 
  sample_info %>% 
  dplyr::filter(Mutation == "Ctrl", Treatment == "UT") %>% 
  .$Sample

depreciation_curve_plot <-
  transcripts_depreciation_curve(sqanti_output = sqanti_class, 
                                 gene_name = Gene,
                                 samples = samples_to_plot)

# Save data ---------------------------------------------------------------
ggsave(plot = depreciation_curve_plot,
       filename = paste0(Gene, "_depreciation_curve_plot_ctrl_UT.png"), 
       path = here::here("results", Gene), 
       width = 6, 
       height = 4, 
       dpi = 600
)

ggsave(plot = depreciation_curve_plot,
       filename = paste0(Gene, "_depreciation_curve_plot_ctrl_UT.svg"), 
       path = here::here("results", Gene), 
       width = 6, 
       height = 4, 
       dpi = 600
)
