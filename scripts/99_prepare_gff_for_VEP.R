# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(rtracklayer)
library(plyranges)

# Load data ---------------------------------------------------------------

lr <- rtracklayer::import(here::here("raw_data", "PB_iPSC_04062021", "SNCA_corrected.gtf.cds.gff"))
transcript <- read.table(here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"), header = TRUE, sep="\t")

# Main --------------------------------------------------------------------

transcripts_to_include <- transcript[transcript$coding == "coding", ]

# Only include transcripts in manuscript
lr_filtered <- lr[lr$transcript_id %in% transcripts_to_include$isoform,]

# change mcols to match required input for VEP
lr_final <- lr_filtered %>% 
  plyranges::bind_ranges(., 
                         GRanges(
                           seqnames = unique(seqnames(lr_filtered)),
                           ranges = IRanges(
                             start = min(start(lr_filtered)),
                             end = max(end(lr_filtered))), 
                           strand = unique(strand(lr_filtered)),
                           type = "gene",
                           source = "PacBio")
  ) %>% 
  plyranges::mutate(Parent = case_when(type == "gene" ~ NA_character_,
                                       type == "transcript" ~ 
                                       TRUE ~ NA_character_))
  
  plyranges::mutate(Parent = case_when(
    type == "gene" ~ NA_character_, # No parent for gene
    type == "transcript" ~ paste0("gene", seqnames, "_", row_number()), # Parent for transcript is gene
    type == "exon" ~ paste0("transcript", seqnames, "_", row_number()), # Parent for exon is transcript
    type == "CDS" ~ paste0("transcript", seqnames, "_", row_number()), # Parent for CDS is transcript
    TRUE ~ NA_character_ # In case there are other types
  ))


lr_final <- lr_filtered %>% 
  plyranges::mutate(ID = gene_id,
                    Name = "SNCA",
                    Parent = transcript_id) %>% 
  plyranges::select(source, type, ID, Name, Parent)


# Save data ---------------------------------------------------------------

rtracklayer::export(lr_final, here::here("results", "SNCA_iPSC.gff"), format = "gff3")
