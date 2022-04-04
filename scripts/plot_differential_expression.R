# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggforce)
library(patchwork)

# Load data ---------------------------------------------------------------

sqanti_class <-
  list(
    iPSC = read.table(here::here("docs", "iPSC", "SNCA_iPSC_classification_filtered.txt"), header = TRUE, sep="\t"),
    Brain = read.table(here::here("docs", "Brain_regions", "SNCA_RFC1_classification_filtered.txt"), header = TRUE, sep="\t")
  )

# Main --------------------------------------------------------------------


iPSC <-
  sqanti_class$iPSC %>% 
  rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
              starts_with('NFLR.Clontech_5p..')) %>% 
  
  # remove ASO samples
  dplyr::select(!ends_with("_M")) %>% 
  
  dplyr::select(isoform, 
                Isoform_class, 
                starts_with(c("045", "828", "829", "831", "C1", "C2", "C3", "C4"))) %>% 
  
  pivot_longer(!c(Isoform_class, isoform), names_to = "Sample", values_to = "count") %>%
  
  dplyr::mutate(Genotype = case_when(grepl("045|828|829|831", Sample) ~ "Mutant",
                                     grepl("C1|C2|C3|C4", Sample) ~ "Control"))
  

## Plot per transcript category ##

iPSC_transcript_category_differential_expression <-
  aggregate(count ~ Isoform_class + Sample + Genotype, 
            data = iPSC, 
            FUN = "sum") %>% 
  
  ggplot(aes(x= factor(Isoform_class, 
                       levels = c("Coding known (complete match)",
                                  "Coding known (alternate 3/5 end)",
                                  "Coding novel",
                                  "Non-coding novel")),
             y = count, 
             fill = Genotype)) + 
  geom_boxplot(colour = "black") +
  scale_fill_brewer(palette = "Dark2") +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..))) +
  labs(x = "Transcript category", y = "Relative Expression (%)") +
  theme_bw() +
  scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                              "Coding known (complete match)" = "Coding known\n(complete match)")) +
  theme(legend.title=element_blank(),
        legend.position="top", legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 12, angle = 45, hjust=1),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1))

## Plot per transcript ##

iPSC_transcript_differential_expression <-
  iPSC %>% 
  dplyr::mutate(isoform = ifelse(isoform == "PB.3.56", "ENST00000336904.7", isoform)) %>%
  dplyr::filter(isoform %in% c("ENST00000336904.7", "PB.3.49")) %>% 
  
  ggplot(aes(x= isoform,
             y = count, 
             fill = Genotype)) + 
  geom_boxplot(colour = "black", show.legend = F) +
  scale_fill_brewer(palette = "Dark2") +
  stat_compare_means(paired = F, 
                     aes(label = paste0("p = ", ..p.format..))) +
  labs(x = "Transcript", y = "Relative Expression (%)") +
  theme_bw() +
  facet_wrap(. ~ factor(Isoform_class, levels = c("Coding known (complete match)", 
                                                  "Coding known (alternate 3/5 end)")), scales = "free", nrow = 1) +
  theme(strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 11, face = "bold"))

# Save data ---------------------------------------------------------------

ggsave(
  plot = iPSC_transcript_category_differential_expression / iPSC_transcript_differential_expression, 
  filename = "iPSC_differential_expression_plot.png", 
  path = here::here("results"), 
  width = 7, 
  height = 8, 
  dpi = 600
)
