library(tidyverse)
library(here)
library(readxl)
library(ggtranscript)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)

# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"),
    path_to_gff = here::here("raw_data", "PB_iPSC_04062021", "SNCA_corrected.gtf.cds.gff"),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx"),
    path_to_ASO = here::here("raw_data", "ASO_design.csv")
  )

# Load data ---------------------------------------------------------------

transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

gff <-
  rtracklayer::import(args$path_to_gff)

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# The probe design
ASO_design <-
  read.table(here::here("raw_data", "ASO_design.csv"), header = T, sep = ",") %>% 
  GRanges()

# Functions ---------------------------------------------------------------

plot_ASO_design <- function(transcripts, gff, ASO, start, end) {
  
  # fill colour to use
  fill_colour <- c("Coding known (complete match)" = "#045a8d",
                   "Coding known (alternate 3/5 end)" = "#74add1",
                   "Coding novel" = "#4d9221",
                   "NMD novel" = "#d53e4f",
                   "Non-coding known" = "#b2abd2",
                   "Non-coding novel" = "#d8daeb")
  
  # this is to define xlim so that plots match
  locus_subset <- 
    data.frame(start = start,
               end = end)
  
  # prepare gff object as a data frame and filter to only include transcripts defined in the transcripts object
  gff <- gff[gff$transcript_id %in% transcripts$isoform, ] %>% data.frame() 
    
  
  # prepare ggtranscript input
  exons <- gff[gff$type == "exon", ]
  introns <- exons %>% to_intron(group_var = "transcript_id")
  cds <- gff[gff$type == "CDS", ]
  
  # transcript plot
  transcript_plot <- 
    exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(height = 0.25) +
    geom_range(
      data = cds
    ) +
    geom_intron(
      data = introns,
      arrow.min.intron.length = 3500,
      strand = "-"
    ) +
    labs(
      y = "Transcript ID",
      x = ""
    ) +
    #xlim(locus_subset$start, locus_subset$end) +
    facet_zoom(xlim = c(locus_subset$start, locus_subset$end), zoom.size = 1.4) +
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  # plot probes
  # Create a new variable 'row' to specify the row for each rectangle
  ASO$y <- rank(-ASO$ASO) / max(rank(-ASO$ASO))
  
  ASO_plot <-
    ASO %>%
    data.frame() %>% 
    ggplot(aes(xmin = start, xmax = end, ymin = y - 0.015, ymax = y + 0.015)) +
    geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8, position = position_nudge(x = 0.05)) +
    xlim(locus_subset$start, locus_subset$end) +
    labs(
      y = "ASO",
      x = "Genomic position (hg38)"
    ) +
    geom_text(aes(label = str_c("ASO ", ASO), x = start - 200, y = y)) +
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

    
  final_plot <- plot_grid(transcript_plot,
                          ASO_plot,
                          ncol = 1,
                          align = "hv", 
                          rel_heights = c(1.6, 0.3), 
                          axis = "lr",
                          label_size = 18)
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# transcripts to include
ASO_design_plot <- plot_ASO_design(transcripts = transcripts, 
                                   gff = gff, 
                                   ASO = ASO_design,
                                   start = min(start(ASO_design)) - 1000, 
                                   end = max(end(ASO_design)) + 1000)

# Save data -------------------------------------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = ASO_design_plot,
    filename = paste0("03a_ASO_design_plot.", ext),
    path = here::here("results", "figures"),
    width = 16,
    height = 18,
    dpi = 600
  )
}
