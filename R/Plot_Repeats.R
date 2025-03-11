PlotRepeats <- function(
    data,
    window_size = 200000,  # Size of each bin (bp)
    step_size   = 200000   # Step between bin starts (bp)
) {
  # Get the list of unique chromosomes
  chr_names <- unique(data$Chromosome)
  
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
    repeat_counts <- sapply(windows, function(w) {
      sum(repeat_midpoints >= w & repeat_midpoints < (w + window_size))
    })
    
    # Create a data frame with the binned results
    results_df <- data.frame(
      Window_Start = windows,
      Repeat_Count = repeat_counts,
      Chromosome   = chr
    )
    
    # Build the plot (no rolling average)
   p <- ggplot(results_df, aes(x = Window_Start, y = Repeat_Count)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = Repeat_Count), fill = "lightblue", alpha = 0.4) +
  labs(
    title = paste("Chromosome:", chr),
    x     = "Genomic Position (binned)",
    y     = "Repeat Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )
    
    # Display the plot for the current chromosome
    print(p)
  }
}
