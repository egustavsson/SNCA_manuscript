library(tidyverse)
library(here)
library(ggtranscript)
library(GenomicRanges)
library(rtracklayer)
library(ggforce)


# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_transcripts = here::here("raw_data", "PB_iPSC_04062021", "SNCA_classification_filtered.txt"),
    path_to_gff = here::here("raw_data", "PB_iPSC_04062021", "SNCA_corrected.gtf.cds.gff"),
    path_to_samples = here::here("raw_data", "PB_iPSC_04062021", "James_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

gff <-
  rtracklayer::import(args$path_to_gff)

Samples <-
  read_excel(args$path_to_samples) %>% 
  data.frame() 

# Functions ---------------------------------------------------------------

plot_transcripts <- function(transcripts, gff, samples, genotype, treatment) {
  
  samples_to_exclude <- samples[samples$Treatment != treatment | samples$Mutation != genotype, ]$Sample
  
  transcripts <- 
    transcripts %>% 
    rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                starts_with('NFLR.Clontech_5p..')) %>% 
    dplyr::select(!samples_to_exclude)
  
  
  # fill colour to use
  fill_colour <- c("Coding known (complete match)" = "#045a8d",
                   "Coding known (alternate 3/5 end)" = "#74add1",
                   "Coding novel" = "#4d9221",
                   "NMD novel" = "#d53e4f",
                   "Non-coding known" = "#b2abd2",
                   "Non-coding novel" = "#d8daeb")
  
  # prepare gff object as a data frame and filter to only include transcripts defined in the transcripts object
  gff <- 
    gff[gff$transcript_id %in% transcripts$isoform, ] %>% 
    data.frame() %>% 
    dplyr::left_join(., transcripts[,c("isoform", "Isoform_class")], by = c("transcript_id" = "isoform"))
  
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
    geom_range(aes(fill = Isoform_class),
               height = 0.25) +
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
    scale_fill_manual(values = fill_colour) +
    ggforce::facet_col(vars(Isoform_class), scale = "free_y", space = "free") +
    theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"))
  
  return(transcript_plot)
  
}

# Main --------------------------------------------------------------------

transcript_plot <- plot_transcripts(transcripts = transcripts, 
                                    gff = gff,
                                    samples = Samples, 
                                    genotype = "Ctrl",
                                    treatment = "UT")



# Save data -------------------------------------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = transcript_plot,
    filename = paste0("01c_transcript_structure_plot.", ext),
    path = here::here("results", "figures"),
    width = 14,
    height = 16,
    dpi = 600
  )
}
