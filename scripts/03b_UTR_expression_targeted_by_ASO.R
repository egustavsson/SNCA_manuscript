library(tidyverse)
library(readxl)
library(here)
library(ggtranscript)
library(GenomicRanges)
library(plyranges)
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

transcripts <- read.table(args$path_to_transcripts, header = TRUE, sep="\t")
gff <- rtracklayer::import(args$path_to_gff)
Samples <- read_excel(args$path_to_samples) %>% data.frame() 
ASO_design <- read.table(here::here("raw_data", "ASO_design.csv"), header = T, sep = ",") %>% GRanges()

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

# Function to get 3'UTRs for each transcript in a GRanges object
get_3UTR <- function(gff) {
  
  # split by transcript so that easier to add conditional statements such as CDS etc.
  split_gff <- gff %>% split(., mcols(.)$transcript_id)
  
  utr_list <- lapply(split_gff, function(x) {
    
    sorted_exons <- x[x$type %in% c("exon", "CDS"), ] %>% sort()
    
    if (NROW(sorted_exons) == 1) {
      last_exon <- x[x$type == "exon", ]
    } else {
      last_exon <- GenomicRanges::setdiff(sorted_exons[1], sorted_exons[2:length(sorted_exons)])
      mcols(last_exon) <- mcols(sorted_exons[1])
    }
    
    return(last_exon)
  })
  
  UTR_gff <- bind_ranges(utr_list)
  
  # Get unique ranges from UTR_gff
  unique_ranges <- unique(UTR_gff)
  
  # Add a column named "UTR" and number it from 1
  unique_ranges$UTR <- seq_along(unique_ranges)
  
  # Compare all ranges to unique_ranges and add "UTR" column
  for (i in seq_along(utr_list)) {
    overlap_idx <- match(utr_list[[i]], unique_ranges)
    utr_list[[i]]$UTR <- unique_ranges$UTR[overlap_idx]
  }
  
  UTR_gff <- plyranges::bind_ranges(utr_list)
  
  return(UTR_gff)
}


# Function to plot ASO design and expression
plot_ASO_design <- function(transcripts, gff, ASO, start, end, samples, genotype, treatment) {
  
  # this is to define xlim so that plots match
  locus_subset <- 
    data.frame(start = start,
               end = end)
  
  # keep gff as GRanges and prepare gff object as a data frame and filter to only include transcripts defined in the transcripts object
  gff_gr <- gff
  gff <- data.frame(gff)
  
  # prepare ggtranscript input
  exons <- gff[gff$type == "exon", ]
  
  # transcript plot
  transcript_plot <- 
    exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = factor(UTR, levels = order),
      fill = overlap_status
    )) +
    geom_range(height = 0.50) +
    scale_fill_manual(
      values = c("in annotation" = "turquoise", "not in annotation" = "grey")
    ) +
    coord_cartesian(xlim = c(locus_subset$start, locus_subset$end)) +
    labs(
      title = "Unique 3' UTRs of SNCA detected",  # Updated title
      y = "Unique 3' UTR",
      x = "Genomic position (hg38)"
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.title=element_blank(),
      legend.position = c(0.5, 0.9),  # Position the legend at the top
      legend.box = "horizontal",  # Display legend in one line
      legend.box.background = element_rect(color = "black", size = 1),
      legend.text = element_text(size = 9),  # Adjust legend text size
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center and style title
    )
  
  # plot probes
  # Create a new variable 'row' to specify the row for each rectangle
  ASO$y <- rank(-ASO$ASO) / max(rank(-ASO$ASO))
  
  ASO_plot <-
    ASO %>%
    data.frame() %>% 
    ggplot(aes(xmin = start, xmax = end, ymin = y - 0.025, ymax = y + 0.025)) +
    geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8, position = position_nudge(x = 0.05)) +
    xlim(locus_subset$start, locus_subset$end) +
    labs(
      title = "ASO designs to target SNCA expression",  # Updated title
      y = "ASO",
      x = "Genomic position (hg38)"
    ) +
    geom_text(aes(label = str_c("ASO ", ASO), x = start - 200, y = y)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  # transcript expression 
  # Filter samples to include
  samples_to_include <- samples[samples$Treatment == treatment & samples$Mutation %in% genotype, ]$Sample
  
  # merge transcripts with above information
  transcripts_UTR <-
    transcripts %>% 
    dplyr::left_join(., data.frame(mcols(UTR_gff)[,c("transcript_id", "UTR")]), by = c("isoform" = "transcript_id"))
  
  
  transcript_UTR_count <-
    transcripts_UTR %>% 
    rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                starts_with('NFLR.Clontech_5p..')) %>% 
    dplyr::select(UTR, all_of(samples_to_include)) %>%
    pivot_longer(!c(UTR), 
                 names_to = "Sample", 
                 values_to = "count") %>% 
    aggregate(data = ., count ~ UTR + Sample, FUN = sum) %>% 
    dplyr::left_join(., 
                     UTR_final %>% data.frame() %>% dplyr::select(UTR, overlap_status), by = "UTR")
  
  # order to plot
  order <- transcript_UTR_count %>% 
    aggregate(data = ., count ~ UTR, FUN = mean) %>% 
    arrange(desc(count)) %>% .$UTR
  
  # transcript expression plot
  expression_plot <-
    transcript_UTR_count %>% 
    ggplot(aes(x = count, y = factor(UTR, levels = order))) +
    geom_boxplot(outlier.shape = NA, aes(fill = overlap_status))+
    geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.25) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      title = "Expression per UTR",  # Updated title
      x = "Relative expression (%)",
      y = ""
    ) +
    scale_fill_manual(
      values = c("in annotation" = "turquoise", "not in annotation" = "grey")
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
  
  # Define what UTRs are targeted by each ASO
  # First need to know which UTR each ASO is targeting
  overlap_results <- findOverlaps(ASO, gff_gr, ignore.strand = TRUE) %>% data.frame()
  
  # loop to get expression targeted by ASO per sample
  expression_by_ASO <- data.frame()  # Initialize an empty data frame to store the results
  
  for (i in unique(overlap_results$queryHits)) {
    
    UTRs <- overlap_results[overlap_results$queryHits == i,]$subjectHits
    
    expression_by_i <- 
      transcript_UTR_count[transcript_UTR_count$UTR %in% UTRs, ] %>% 
      aggregate(data = ., count ~ Sample, FUN = sum) %>% 
      dplyr::mutate(ASO = i)
    
    # Append the results to the data frame
    expression_by_ASO <- rbind(expression_by_ASO, expression_by_i)
  }
  
  # Percentage expression targeted by ASO plot
  expression_by_ASO_plot <-
    expression_by_ASO %>% 
    ggplot(aes(x = factor(ASO),
               y = count)) +
    geom_boxplot(outlier.shape = NA, fill = "grey") +
    geom_jitter(color="black", size=0.5, alpha=0.9, width = 0.25) +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       limits = c(0, 100),
                       expand = c(0, 0)) +
    labs(
      title = "SNCA expression targeted by each ASO",  # Updated title
      x = "ASO",
      y = "Relative expression (%)") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size=16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"))
  
  
  # Final plot
  final_plot <- plot_grid(
    transcript_plot, 
    expression_plot, 
    ASO_plot, 
    expression_by_ASO_plot, 
    ncol = 2, 
    nrow = 2, 
    align = "hv", 
    axis = "tblr", 
    rel_heights = c(1.6, 1), labels = c("A", "B", "C", "D"), 
    label_size = 18
  )
  
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# Get 3' UTRs from long read data
UTR_gff <- gff[gff$transcript_id %in% transcripts$isoform, ] %>% get_3UTR()

# Get SNCA reference 3' UTRs
ref_UTR_gff <- ref[ref$gene_name == "SNCA"] %>% get_3UTR()


overlaps <- findOverlaps(UTR_gff, ref_UTR_gff, ignore.strand = TRUE)

# Initialize a new column in UTR_gff
UTR_gff$overlap_status <- "not in annotation"

# Mark the ranges that overlap as "in annotation"
UTR_gff$overlap_status[subjectHits(overlaps)] <- "in annotation"


UTR_final <-
  UTR_gff %>% 
  plyranges::select(type,
                    UTR,
                    overlap_status) %>% unique()

# transcripts to include
ASO_design_plot <- plot_ASO_design(transcripts = transcripts, 
                                   gff = UTR_final, 
                                   ASO = ASO_design,
                                   samples = Samples, 
                                   genotype = "Ctrl", 
                                   treatment = "UT",
                                   start = min(start(ASO_design)) - 1000, 
                                   end = max(end(ASO_design)) + 500)

# Save data -------------------------------------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = ASO_design_plot,
    filename = paste0("03b_UTR_expression_ASO_design_plot.", ext),
    path = here::here("results", "figures"),
    width = 16,
    height = 12,
    dpi = 600, bg = "white"
  )
}
