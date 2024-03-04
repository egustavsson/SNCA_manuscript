# This script will aggregate transcripts expression by calculating a mean across indictions

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(stringr)

# Load data ---------------------------------------------------------------

transcript_expression_data <- 
  read.table(here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"), header = TRUE, sep="\t")

sample_data <- 
  read_excel(here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")) %>% 
  data.frame()

# Main --------------------------------------------------------------------

aggregated_mean_transcript_expression <-
  transcript_expression_data %>% 
  rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
            starts_with('NFLR.Clontech_5p..')) %>% 
  dplyr::select(Isoform_class, isoform, all_of(sample_data$Sample)) %>%
  pivot_longer(!c(Isoform_class, isoform), 
               names_to = "Sample", 
               values_to = "count") %>% 
  dplyr::left_join(., 
                   dplyr::select(sample_data, 
                                 Sample, 
                                 Line,
                                 Mutation, 
                                 Treatment), 
                   by = c("Sample" = "Sample")) %>% 
  aggregate(count ~ Isoform_class + isoform + Mutation + Line + Treatment,
            data = .,
            FUN = "mean") %>% 
  dplyr::mutate(case_ctrl = ifelse(Mutation == "Ctrl", "Control", "Mutant")) 

