# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggforce)

# Load data ---------------------------------------------------------------

sqanti_class <-
  list(
    iPSC = read.table(here::here("docs", "iPSC", "SNCA_iPSC_classification_processed.txt"), header = TRUE, sep="\t"),
    Brain = read.table(here::here("docs", "Brain_regions", "SNCA_RFC1_classification_processed.txt"), header = TRUE, sep="\t")
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
  
  dplyr::mutate(Genotype = case_when(grepl("045", Sample) ~ "G51D",
                                     grepl("828|829", Sample) ~ "A53T",
                                     grepl("831", Sample) ~ "SNCAx3",
                                     grepl("C1|C2|C3|C4", Sample) ~ "Control"))



## Unique transcripts
iPSC_transcripts_plot <-
  iPSC %>% 
  select(!c(Sample, count, Genotype)) %>%
  unique() %>%
  dplyr::mutate(Tissue = "iPSC mDA neurons") %>% 
  
  ggplot(aes(x = factor(Isoform_class, 
                        levels = c("Coding known (complete match)",
                                   "Coding known (alternate 3/5 end)",
                                   "Coding novel",
                                   "Non-coding novel"))
             , fill = Isoform_class)) +
  geom_bar(show.legend = F, colour = "Black") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                              "Coding known (complete match)" = "Coding known\n(complete match)")) +
  labs(y = "No. unique transcripts", x = "Transcript category") +
  facet_wrap(~ Tissue) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, 
                                   hjust=1,
                                   size = 10))



## Unique transcripts per genotype ##  
iPSC_transcripts_per_genotype_plot <-
  
  iPSC %>% 
  select(!c(Sample, count)) %>%
  unique() %>% 
  
  ggplot(aes(x = factor(Isoform_class, 
                        levels = c("Coding known (complete match)",
                                   "Coding known (alternate 3/5 end)",
                                   "Coding novel",
                                   "Non-coding novel"))
             , fill = Isoform_class)) +
  geom_bar(show.legend = F, colour = "Black") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                              "Coding known (complete match)" = "Coding known\n(complete match)")) +
  labs(y = "No. unique transcripts", x = "Transcript category") +
  facet_wrap(~ factor(Genotype, levels = c("Control",
                                           "A53T",
                                           "G51D",
                                           "SNCAx3"))) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, 
                                   hjust=1,
                                   size = 10))

# Save data ---------------------------------------------------------------

ggsave(
  plot = iPSC_transcripts_plot, 
  filename = "iPSC_unique_transcripts_plot.png", 
  path = here::here("results"), 
  width = 5, 
  height = 5, 
  dpi = 600
)

ggsave(
  plot = iPSC_transcripts_per_genotype_plot, 
  filename = "iPSC_unique_transcripts_per_genotype_plot.png", 
  path = here::here("results"), 
  width = 6, 
  height = 5, 
  dpi = 600
)
