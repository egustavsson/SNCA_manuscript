# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(patchwork)
library(rtracklayer)
library(plyranges)
library(here)

# Load data ---------------------------------------------------------------

Filtered_tx <-
  list(
    iPSC = read.table(here::here("docs", "iPSC", "SNCA_iPSC_classification_filtered.txt"), header = TRUE, sep="\t"),
    Brain = read.table(here::here("docs", "Brain_regions", "SNCA_RFC1_classification_filtered.txt"), header = TRUE, sep="\t")
  )

LR_gff <- 
  list(
    iPSC = rtracklayer::import(here::here("docs", "iPSC", "SNCA_iPSC_corrected.gtf.cds.gff")),
    Brain = rtracklayer::import(here::here("docs", "Brain_regions", "SNCA_RFC1_corrected.gtf.cds.gff"))
  )

# Functions ---------------------------------------------------------------

get_TOI <- 
  
  function(gff, pb_ids) {
    
    gff <- gff[gff$transcript_id %in% pb_ids]
    gff_exons <- gff[gff$type == "exon"]
    gff_cds <- gff[gff$type == "CDS"]
    
    # return as list as we need both exons and cds
    gff_exons_cds <- 
      list(
        exons = data.frame(gff_exons),
        cds = data.frame(gff_cds)
      )
    
    return(gff_exons_cds)
  }


get_MANE <- 
  
  function(ref, MANE_id) {
    
    # remove any NA transcript ids (i.e. type == "gene")
    MANE <- ref[!is.na(ref$transcript_id)]
    
    # remove to .XX after ENST
    GenomicRanges::mcols(MANE)[["transcript_id"]] <- 
      GenomicRanges::mcols(MANE)[["transcript_id"]] %>%
      stringr::str_remove("\\..*")
    
    MANE <- MANE[MANE$transcript_id == mane_id, ]
    MANE_exons <- MANE[MANE$type == "exon"]
    MANE_cds <- MANE[MANE$type == "CDS"]
    
    MANE_exons_cds <- 
      list(
        exons = MANE_exons,
        cds = MANE_cds
      )
    
    return(MANE_exons_cds)
  }

annotate_TX <- 
  
  function(tx) {
    
    sapply(tx, function(x){
      tmp <- Filtered_tx$iPSC
      x$transcript_biotype <- tmp[match(x$transcript_id, tmp$isoform),]$Isoform_class
      
      return(x)
    }, simplify=F)
    
  }

# Main --------------------------------------------------------------------

# All filtered SNCA transcripts #

SNCA_iPSC <- 
  get_TOI(
    gff = LR_gff$iPSC,
    pb_ids = Filtered_tx$iPSC$isoform) %>% 
  annotate_TX(tx = .) %>% 
  
  ## this needs to be added as a conditional statement for more general use ##
  lapply(., function(x){
    x %>% 
      dplyr::mutate(transcript_id = ifelse(transcript_id == "PB.3.56", "ENST00000336904.7", transcript_id)) })
  

  
  
SNCA_exons <- 
  SNCA_iPSC$exons 
  

SNCA_introns <- 
  SNCA_exons %>% 
  to_intron(group_var = "transcript_id")

SNCA_cds <- 
  SNCA_iPSC$cds

SNCA_base_plot <-
  SNCA_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = transcript_biotype),
    height = 0.25 
  ) +
  geom_range(
    data = SNCA_cds
  ) +
  geom_intron(
    data = SNCA_introns,
    arrow.min.intron.length = 3500,
    strand = "-"
  ) +
  labs(
    y = "Transcript name",
    x = "Genomic position"
  ) +
  scale_fill_brewer(palette = "Dark2")

SNCA_base_plot

## SNCA all filtered rescaled ##

SNCA_rescaled <-
  shorten_gaps(SNCA_exons, SNCA_introns, group_var = "transcript_id")

SNCA_rescaled_base_plot <-
  SNCA_rescaled %>% 
  dplyr::filter(type == "exon") %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(aes(fill = transcript_biotype)) +
  geom_intron(
    data = SNCA_rescaled %>% dplyr::filter(type == "intron"),
    arrow.min.intron.length = 300,
    strand = "-"
  ) +
  labs(
    y = "Transcript name",
    x = "Rescaled position"
  ) +
  scale_fill_discrete(name = "Transcript biotype") +
  scale_fill_brewer(palette = "Dark2") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

legend <- ggpubr::as_ggplot(ggpubr::get_legend(SNCA_rescaled_base_plot))

SNCA_base <- ggpubr::ggarrange(legend, 
                               SNCA_base_plot + theme(legend.position = "none") + SNCA_rescaled_base_plot + theme(legend.position = "none",
                                                                                                                  axis.title.y = element_blank(),
                                                                                                                  axis.text.y = element_blank(),
                                                                                                                  axis.ticks.y = element_blank()), 
                               nrow = 2,
                               heights = c(0.1, 0.9))

# Save data ---------------------------------------------------------------

ggsave(
  plot = SNCA_base, 
  filename = "SNCA_base_plot.png", 
  path = here::here("results"), 
  width = 12, 
  height = 7, 
  dpi = 600, 
  bg = "white"
)

