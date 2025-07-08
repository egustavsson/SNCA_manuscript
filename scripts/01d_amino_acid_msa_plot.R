# Load libraries ----------------------------------------------------------

library(msa)
library(Biostrings)
library(here)
library(dplyr)
library(magrittr)

# Load data ---------------------------------------------------------------

iPSC <-
  read.table(here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"), header = TRUE, sep="\t")


# Main --------------------------------------------------------------------

seqs <- iPSC %>% 
  dplyr::select(ORF_length, ORF_seq) %>% 
  na.omit() %>% 
  unique() %>% 
  dplyr::mutate(ORF_length = paste0(ORF_length, "aa")) %>% 
  arrange(desc(ORF_length))

seq_vec <- setNames(seqs$ORF_seq, seqs$ORF_length)


seqs_AA <- AAStringSet(seq_vec)

# Perform multiple sequence alignment
alignment <- msa(seqs_AA, method = "ClustalOmega", order = "input", type = "protein")  

# Optional: create a PDF or image
msaPrettyPrint(alignment,
               file = "SNCA_isoforms_alignment.pdf",
               output = "pdf",
               showNames = "left",
               showLogo = "none",
               shadingMode="identical",
               askForOverwrite = FALSE,
               verbose = FALSE, 
               showConsensus = "none", 
               showNumbering = "none",
               showLegend = FALSE, 
               paperWidth = 15, 
               paperHeight = 6, 
               furtherCode=c(
                 "\\defconsensus{.}{lower}{upper}",
                 "\\showruler{1}{top}")
              )
