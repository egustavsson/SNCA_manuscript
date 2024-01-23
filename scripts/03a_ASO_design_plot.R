library(tidyverse)
library(here)
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
  
  # this is to define xlim so that plots match
  locus_subset <- 
    data.frame(start = start,
               end = end)
  
  # prepare ASO object as a data frame
  ASO_design <- data.frame(ASO_design)
  
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
    xlim(locus_subset$start, locus_subset$end) +
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  # plot probes
  # Create a new variable 'row' to specify the row for each rectangle
  ASO$y <- runif(nrow(data.frame(ASO_design)))
  
  
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
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  final_plot <- plot_grid(transcript_plot,
                          ASO_plot,
                          ncol = 1,
                          align = "hv", 
                          rel_heights = c(1, 0.5), 
                          axis = "lr",
                          label_size = 18)
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# transcripts to include
transcript_plot <- plot_ASO_design(transcripts = transcripts, 
                                   gff = gff, 
                                   ASO = ASO_design,
                                   start = 89699710, 
                                   end = 89838977)

transcript_plot_zoomed <- plot_ASO_design(transcripts = transcripts, 
                                          gff = gff, 
                                          ASO = ASO_design,
                                          start = min(start(ASO_design)) - 1000, 
                                          end = max(end(ASO_design)) + 1000)


plot_panel <- 
  plot_grid(transcript_plot,
            transcript_plot_zoomed, 
            nrow = 1,
            align = 'hv',
            labels = c('A', 'B'))
  
# Save data -------------------------------------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = plot_panel,
    filename = paste0("03a_ASO_design_plot.", ext),
    path = here::here("results", "figures"),
    width = 18,
    height = 16,
    dpi = 600
  )
}
