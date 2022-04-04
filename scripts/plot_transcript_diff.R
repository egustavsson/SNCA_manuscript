# Load libraries ----------------------------------------------------------

library(ggtranscript)
library(tidyverse)
library(here)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)

# Load data ---------------------------------------------------------------

TOI <- 
  read.table(
    here::here("docs", "iPSC", "SNCA_iPSC_classification_filtered.txt"), 
    header = TRUE, 
    sep="\t")
 
lr <- 
  rtracklayer::import(
    here::here("docs", "iPSC", "SNCA_iPSC_corrected.gtf.cds.gff")
    ) %>%
  # Only include filtered transcripts
  plyranges::filter(transcript_id %in% TOI$isoform)

# also download and load reference annotation 
# ref_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")
# 
# if(!file.exists(ref_path)) {
#   
#   download.file(
#     url = paste0(
#       "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
#       "gencode.v38.annotation.gtf.gz"
#     ),
#     destfile = ref_path
#   )
#   
# }
# 
# R.utils::gunzip(ref_path, remove = TRUE)
# 
# ref <- rtracklayer::import(stringr::str_remove(ref_path, "\\.gz"))

ref <- 
  rtracklayer::import("/home/egust/Projects/pseudogenes/raw_data/GENCODE_annotations/gencode.v38.annotation.gtf.gz")

# Functions ---------------------------------------------------------------

get_lr_tx_of_interest <- 
  
  function(lr, pb_ids) {
    
    lr <- lr[lr$transcript_id %in% pb_ids]
    lr_exons <- lr[lr$type == "exon"]
    lr_cds <- lr[lr$type == "CDS"]
    
    # return as list as we need both exons and cds
    lr_exons_cds <- 
      list(
        exons = lr_exons,
        cds = lr_cds
        )
    
    return(lr_exons_cds)
  }

get_mane <- 
  
  function(ref, mane_id) {
    
    # remove any NA transcript ids (i.e. type == "gene")
    mane <- ref[!is.na(ref$transcript_id)]
    
    # remove to .XX after ENST
    GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>%
      stringr::str_remove("\\..*")
    
    mane <- mane[mane$transcript_id == mane_id, ]
    mane_exons <- mane[mane$type == "exon"]
    mane_cds <- mane[mane$type == "CDS"]
    
    mane_exons_cds <- 
      list(
        exons = mane_exons,
        cds = mane_cds
        )
    
    return(mane_exons_cds)
  }

plot_diff <- function(lr_exons_cds, 
                      mane_exons_cds, 
                      lr_mane_diffs, 
                      MANE_canonical = "MANE"
) {
  
  # merge mane and lr data and convert to data.frame() for plotting
  # convert transcript_id to factor to make sure mane is at top
  transcript_order <- c(
    lr_exons_cds$exons$transcript_id %>% unique(),
    mane_exons_cds$exons$transcript_id %>% unique()
  )
  
  lr_mane_exons_df <- c(lr_exons_cds$exons, mane_exons_cds$exons) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  lr_mane_cds_df <- c(lr_exons_cds$cds, mane_exons_cds$cds) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  
  # plot diff plot
  diff_plot <- lr_mane_exons_df %>% 
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = transcript_id
    )) + 
    geom_range(
      height = 0.25, 
      fill = "white"
    ) + 
    geom_range(
      data = lr_mane_cds_df
    ) + 
    geom_intron(
      data = to_intron(lr_mane_exons_df, "transcript_id"), 
      aes(strand = strand), 
      arrow.min.intron.length = 400
    ) + 
    geom_range(
      data = lr_mane_diffs, 
      aes(
        fill = diff_type, 
        colour = diff_type
      ), 
      alpha = 0.2, 
      linetype = 2
    ) + 
    scale_y_discrete(name = "Transcript ID") + 
    scale_x_continuous(name = "Genomic position") + 
    scale_fill_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    scale_colour_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    theme_bw() + 
    theme(legend.position = "top")
  
  return(diff_plot)
  
}

# Main --------------------------------------------------------------------

##### SNCA #####

SNCA_lr_exons_cds <- get_lr_tx_of_interest(
  lr = lr, 
  pb_ids = TOI$isoform)

SNCA_mane_exons_cds <- get_mane(ref, mane_id = "ENST00000394991")

# obtain differences between MANE and lr exons
SNCA_lr_mane_diffs <- 
  ggtranscript::to_diff(
    exons = SNCA_lr_exons_cds$exons %>% as.data.frame(), 
    ref_exons = SNCA_mane_exons_cds$exons %>% as.data.frame(), 
    group_var = "transcript_id"
  )

SNCA_lr_mane_diff_plot <- 
  plot_diff(
    lr_exons_cds = SNCA_lr_exons_cds, 
    mane_exons_cds = SNCA_mane_exons_cds, 
    lr_mane_diffs = SNCA_lr_mane_diffs
  )


# Save data ---------------------------------------------------------------

ggsave(
  plot = SNCA_lr_mane_diff_plot, 
  filename = "SNCA_lr_mane_diff_plot.png", 
  path = here::here("results"), 
  width = 6, 
  height = 4, 
  dpi = 600
)

