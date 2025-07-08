# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggtranscript)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(cowplot)

# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_transcripts = here::here("raw_data", "ONT_mouse", "Snca_transcripts_controls_ref_guided_corrected.gtf.cds.gff"),
    path_to_expression = here::here("raw_data", "ONT_mouse", "Transcripts_Snca.csv")
  )

# Load data ---------------------------------------------------------------

transcripts <-
  rtracklayer::import(args$path_to_transcripts)

expression <-
  read.table(args$path_to_expression, header = TRUE, sep=",")

# Functions ---------------------------------------------------------------

# Function to plot transcripts
plot_transcript_scaled <- function(gff) {
  
  # Ensure gff is a data frame
  if (!is.data.frame(gff)) {
    gff <- data.frame(gff)
  }
  
  # Extract exons and CDS
  exons <- gff %>% filter(type == "exon")
  cds <- gff %>% filter(type == "CDS")
  
  # Add UTR regions
  utr_data <- add_utr(exons = exons, cds = cds, group_var = "transcript_id")
  
  # Rescale exons and introns with shorten_gaps
  rescaled_data <- shorten_gaps(
    exons = utr_data,
    introns = to_intron(utr_data, "transcript_id"),
    group_var = "transcript_id"
  )
  
  # Plot
  transcript_plot <- 
    rescaled_data %>%
    dplyr::filter(type == "CDS") %>% 
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(fill = "#5A5A5A") +
    geom_range(
      data = rescaled_data %>% dplyr::filter(type == "UTR"),
      height = 0.25,
      fill = "#F5F5F5"
    ) +
    geom_intron(
      data = to_intron(
        rescaled_data %>% dplyr::filter(type != "intron"),
        "transcript_id"),
      arrow.min.intron.length = 100,
      strand = "-"
    ) +
    labs(
      title = "Snca transcripts",
      y = "Transcript ID",
      x = "Genomic Position (Rescaled)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size=16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold")
    )
  
  return(transcript_plot)
}

# Function to get 3'UTRs for each transcript in a GRanges object
get_3UTR <- function(gff) {
  
  # Split the GFF by transcript_id for easier conditional processing
  split_gff <- gff %>%
    split(mcols(.)$transcript_id)
  
  # Process each transcript to extract the 3' UTR
  utr_list <- lapply(split_gff, function(transcript_ranges) {
    # Filter and sort exons and CDS ranges
    sorted_exons <- transcript_ranges[transcript_ranges$type %in% c("exon", "CDS"), ] %>% 
      sort()
    
    # Handle cases where there is only one range
    if (NROW(sorted_exons) == 1) {
      last_exon <- transcript_ranges[transcript_ranges$type == "exon", ]
    } else {
      # Subtract all other ranges from the first range to isolate the 3' UTR
      last_exon <- GenomicRanges::setdiff(
        sorted_exons[1], 
        sorted_exons[2:length(sorted_exons)]
      )
      # Preserve metadata from the first range
      mcols(last_exon) <- mcols(sorted_exons[1])
    }
    
    return(last_exon)
  })
  
  # Combine the processed 3' UTRs into a single range
  UTR_gff <- bind_ranges(utr_list)
  
  # Extract unique ranges
  unique_ranges <- UTR_gff
  #unique(UTR_gff)
  
  # Find overlaps between unique ranges and the original UTR ranges
  UTR_hits <- findOverlaps(unique_ranges, UTR_gff, type = "equal")
  
  # Annotate unique ranges with transcript IDs
  unique_ranges <- unique_ranges %>%
    plyranges::mutate(
      transcript_sharing_UTR = sapply(
        seq_along(.),
        function(i) {
          # Find indices of overlapping ranges
          hits <- subjectHits(UTR_hits[queryHits(UTR_hits) == i])
          # Combine transcript IDs from overlapping ranges
          paste(unique(mcols(UTR_gff)$transcript_id[hits]), collapse = ",")
        }
      )
    )
  
  return(unique_ranges)
}


# Plot UTRs
plot_UTRs_with_expression <- function(gff, expression) {
  
  # Extract exons and CDS
  exons <- gff %>% 
    unique() %>% # remove UTRs that are identical
    data.frame() %>% 
    filter(type == "exon")
  
  exons_rescaled <- shorten_gaps(
    exons = exons, 
    introns = to_intron(exons, "transcript_id"), 
    group_var = "transcript_id"
  )
  
  
  # transcript plot
  UTR_plot <- 
    exons_rescaled %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_sharing_UTR
    )) +
    geom_range(height = 0.50, fill = "#F5F5F5") +
    labs(
      title = "3' UTRs of Snca",
      y = "Unique 3' UTRs of Snca Transcripts",
      x = ""
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.title=element_blank(),
      legend.position = c(0.5, 0.9),  # Position the legend at the top
      legend.box = "horizontal",  # Display legend in one line
      legend.box.background = element_rect(color = "black", size = 1),
      legend.text = element_text(size = 9),  # Adjust legend text size
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center and style title
    )
  
  # Process expression
  # Normalize columns with "Barcode" in their name
  barcode_columns <- grep("Barcode", names(expression), value = TRUE) 
  expression[barcode_columns] <- lapply(expression[barcode_columns], function(column) {
    column / sum(column) * 100 # 
  })
  
  expression_final <-
    expression %>% 
    dplyr::full_join(., 
                     dplyr::select(
                       data.frame(gff), 
                       transcript_id, 
                       transcript_sharing_UTR),
                     by = "transcript_id") %>% 
    dplyr::select(-transcript_id) %>% 
    pivot_longer(!c(transcript_sharing_UTR), 
                 names_to = "Sample", 
                 values_to = "count") %>%
    aggregate(data = ., count ~ transcript_sharing_UTR + Sample, FUN = sum)
  
  # Plot expression
  expression_plot <-
    expression_final %>% 
    ggplot(aes(x = count, y = transcript_sharing_UTR)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.25) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      title = "Expression per UTR",  # Updated title
      x = "Relative expression (%)",
      y = ""
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(legend.title=element_blank(),
          legend.position = c(0.5, 0.9),  # Position the legend at the top
          legend.box = "horizontal",  # Display legend in one line
          legend.box.background = element_rect(color = "black", size = 1),
          legend.text = element_text(size = 9),
          plot.title = element_text(size=16, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  
  final_plot <- plot_grid(
    UTR_plot,
    expression_plot,
    nrow = 1, 
    align = "h", 
    rel_widths = c(2, 1), 
    label_size = 18
  )
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# Plot mouse Snca transcripts
transcript_scaled <- plot_transcript_scaled(gff = transcripts)

# Get 3' UTRs
mouse_Snca_UTR <- get_3UTR(gff = transcripts)

# Plot UTRs and expression
mouse_UTRs <- plot_UTRs_with_expression(gff = mouse_Snca_UTR, expression = expression)

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_scaled,
  filename = "06c_Mouse_Snca_transcript_plot.svg",
  path = here::here("results", "figures"),
  width = 9,
  height = 6,
  dpi = 600
)

ggsave(
  plot = mouse_UTRs,
  filename = "06c_mouse_UTRs_plot.svg",
  path = here::here("results", "figures"),
  width = 12,
  height = 6,
  dpi = 600
)
