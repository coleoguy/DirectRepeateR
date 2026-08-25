PlotRepeats <- function(
    data,
    window_size = 200000,  # Size of each bin (bp)
    step_size   = 200000   # Step between bin starts (bp)
) {
  # Get the list of unique chromosomes
  chr_names <- unique(data$Chromosome)
  
  plots <- vector("list", length(chr_names))
  names(plots) <- chr_names
  
  # Loop through each chromosome
  for (chr in chr_names) {
    # Subset for the current chromosome
    chr_data <- subset(data, Chromosome == chr)
    
    # Determine the chromosome length
    chr_length <- max(c(chr_data$End, chr_data$Match_End))
    
    # Calculate midpoint of each repeat
    repeat_midpoints <- ((chr_data$Start + chr_data$End) / 2 +
                           (chr_data$Match_Start + chr_data$Match_End) / 2) / 2
    
    # Define bin starts
    windows <- seq(1, chr_length, by = step_size)
    
    # Count repeats in each bin
    if (window_size == step_size) {
      # Non-overlapping bins: a single binning pass counts every midpoint
      repeat_counts <- tabulate(
        findInterval(repeat_midpoints, windows),
        nbins = length(windows)
      )
    } else {
      # Overlapping or gapped bins: count midpoints in [w, w + window_size)
      # via two sorted lookups instead of one pass per window
      sorted_mid <- sort(repeat_midpoints)
      hi <- findInterval(windows + window_size, sorted_mid,
                         left.open = TRUE)
      lo <- findInterval(windows, sorted_mid, left.open = TRUE)
      repeat_counts <- hi - lo
    }
    
    # Create a data frame with the binned results
    results_df <- data.frame(
      Window_Start = windows,
      Repeat_Count = repeat_counts,
      Chromosome   = chr
    )
    
    # Build the plot (no rolling average)
  p <- ggplot(results_df, aes(x = Window_Start, y = Repeat_Count)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = 0, ymax = Repeat_Count), fill = "lightblue", alpha = 0.4) +
  labs(
    title = paste("Chromosome:", chr),
    x     = "Genomic Position (binned)",
    y     = "Repeat Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.line          = element_line(color = "black"),
    axis.ticks         = element_line(color = "black"),
    axis.ticks.length  = grid::unit(0.2, "cm")
  ) 
    # Display the plot for the current chromosome
    print(p)
    plots[[chr]] <- p
  }
  
  # Return the plots (invisibly) so users can save or modify them
  invisible(plots)
}
