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
  
  dplyr::mutate(Genotype = case_when(grepl("045", Sample) ~ "G51D",
                                     grepl("828|829", Sample) ~ "A53T",
                                     grepl("831", Sample) ~ "SNCAx3",
                                     grepl("C1|C2|C3|C4", Sample) ~ "Control"))


## Plot transcript category expression ##

iPSC_transcript_category_expression_plot <-

aggregate(count ~ Isoform_class + Genotype, 
          data = iPSC, 
          FUN = "sum") %>% 
  
  ggplot(aes(x=factor(Genotype, levels = c("SNCAx3",
                                           "G51D",
                                           "A53T",
                                           "Control")), 
             y = count, 
             fill = Isoform_class)) + 
  geom_col(position="fill", colour = "black", show.legend = F) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip() +
  labs(y = "Transcript category expression", x = "Genotype") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

## Plot unique transcript expression ##

iPSC_unique_transcript_expression_plot <-
  iPSC %>%  
  
  group_by(isoform) %>% 
  ggplot(aes(x = reorder(isoform, -count),
             y = count,
             fill = Isoform_class)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "Transcript expression (NE)", x = "Transcripts ranked by expression (NE)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save data ---------------------------------------------------------------

ggsave(
  plot = iPSC_unique_transcript_expression_plot / iPSC_transcript_category_expression_plot, 
  filename = "iPSC_transcript_expression_plot.png", 
  path = here::here("results"), 
  width = 9, 
  height = 6, 
  dpi = 600
)
