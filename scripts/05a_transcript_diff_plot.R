
# Load libraries ----------------------------------------------------------

# as ggtranscript is in development during the creation of this script
# i have included the code to install the exact version using (0.99.2)
# based on the commit SHA
# devtools::install_github("dzhang32/ggtranscript@88d23b43b42de3e49a67dc8130f5263b6dcf81d1")
library(ggtranscript)
library(tidyverse)
library(here)
library(R.utils)
library(rtracklayer)
library(GenomicRanges)
library(ggforce)
library(cowplot)

# Load data ---------------------------------------------------------------

lr <- rtracklayer::import(here::here("raw_data", "PB_iPSC_04062021", "SNCA_corrected.gtf.cds.gff"))
transcript <- read.table(here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"), header = TRUE, sep="\t")

# also download and load reference annotation 
ref_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")

if(!file.exists(ref_path)) {
  
  download.file(
    url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
      "gencode.v38.annotation.gtf.gz"
    ),
    destfile = ref_path
  )
  
}

R.utils::gunzip(ref_path, remove = TRUE)

ref <- rtracklayer::import(stringr::str_remove(ref_path, "\\.gz"))

# Functions ---------------------------------------------------------------

get_lr_tx_of_interest <- function(lr, novel_ids, pb_ids) {
  
  lr <- lr[lr$transcript_id %in% pb_ids]
  # added this so that the names include aa length
  for (i in 1:length(lr)) {
    # Get the transcript_id of the current row
    current_transcript_id <- lr$transcript_id[i]
    
    # Find the corresponding new_id from the replacement table
    new_id <- novel_ids$new_id[novel_ids$isoform == current_transcript_id]
    
    # Check if a matching new_id is found
    if (length(new_id) > 0) {
      # Replace the transcript_id with the new_id
      lr$transcript_id[i] <- new_id
    }
  }
  
  
  lr_exons <- lr[lr$type == "exon"]
  lr_cds <- lr[lr$type == "CDS"]
  
  # return as list as we need both exons and cds
  lr_exons_cds <- list(
    exons = lr_exons, 
    cds = lr_cds
  )
  
  return(lr_exons_cds)
  
}

get_mane <- function(ref, mane_id) {
  
  # remove any NA transcript ids (i.e. type == "gene")
  mane <- ref[!is.na(ref$transcript_id)] 
  
  # remove to .XX after ENST
  GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>% 
    stringr::str_remove("\\..*")
  
  mane <- mane[mane$transcript_id == mane_id, ]
  mane_exons <- mane[mane$type == "exon"]
  mane_cds <- mane[mane$type == "CDS"]
  
  mane_exons_cds <- list(
    exons = mane_exons, 
    cds = mane_cds
  )
  
  return(mane_exons_cds)
  
}

plot_diff <- function(lr_exons_cds, 
                      mane_exons_cds, 
                      lr_mane_diffs, 
                      MANE_canonical = "MANE",
                      zoom_start, 
                      zoom_end
) {
  
  # merge mane and lr data and convert to data.frame() for plotting
  # convert transcript_id to factor to make sure mane is at top
  transcript_order <- c(
    lr_exons_cds$exons$transcript_id %>% unique() %>% sort(),
    mane_exons_cds$exons$transcript_id %>% unique() %>% sort()
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
    scale_x_continuous(name = "Genomic position (hg38)") + 
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
    theme(legend.position = "top",
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 12, face = "bold")
    )
  
  
  # add zoomed plots
  zoom_3end <- diff_plot + coord_cartesian(xlim = c(89724000, 89730000))
  zoom_5end <- diff_plot + coord_cartesian(xlim = c(89821400, 89837400))
  
  final_plot <- 
    plot_grid(
      diff_plot,
      plot_grid(
        zoom_3end + theme(legend.position = "none"),
        zoom_5end + theme(legend.position = "none",
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())),
      nrow = 2, 
      ncol = 1
    )
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# Get all unique ORFs
ORFs <- transcript %>% dplyr::select(isoform, ORF_length) %>% na.omit() %>% dplyr::filter(ORF_length != "140")

# to add aa length into transcripts id i need to generate them to be used in get_lr_tx_of_interest()
ORF_ids <- ORFs %>% dplyr::mutate(new_id = paste0(ORF_length, "aa", " (", isoform, ")")) %>% dplyr::select(isoform, new_id)

SNCA_lr_exons_cds <- get_lr_tx_of_interest(
  lr = lr,
  novel_ids = ORF_ids,
  pb_ids = ORFs$isoform)

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

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = SNCA_lr_mane_diff_plot,
    filename = paste0("05a_SNCA_mane_diff_plot.", ext),
    path = here::here("results", "figures"),
    width = 12,
    height = 8,
    dpi = 600
  )
}


