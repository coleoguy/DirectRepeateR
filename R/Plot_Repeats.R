PlotRepeats <- function(data, window_size = NULL, step_size = NULL) {

  # Set default values if parameters are NULL
  if (is.null(window_size)) {
    window_size <- 20000
  }
  if (is.null(step_size)) {
    step_size <- 20000
  }
  
  # Get unique chromosome names
  chr_names <- unique(data$Chromosome)
  
  # Process each chromosome and generate a plot
  for (chr in chr_names) {
    # Subset data for the current chromosome
    chr_data <- subset(data, Chromosome == chr)
    
    # Get max positional value (as a proxy for chromosome length)
    chr_length <- max(c(chr_data$End_Position, chr_data$Match_End_Position))
    
    # Calculate the midpoint for each repeat
    repeat_midpoints <- ((chr_data$Start_Position + chr_data$End_Position) / 2 +
                           (chr_data$Match_Position + chr_data$Match_End_Position) / 2) / 2
    
    # Create sliding windows
    windows <- seq(1, chr_length, by = step_size)
    
    # Count repeats in each window
    repeat_counts <- sapply(windows, function(w) {
      sum(repeat_midpoints >= w & repeat_midpoints < (w + window_size))
    })
    
    # Store results for the current chromosome
    results_df <- data.frame(
      Window_Start = windows,
      Repeat_Count = repeat_counts,
      Chromosome = chr
    )
    
    # Create the ggplot object
    ggplot_object <- ggplot(results_df, aes(x = Window_Start, y = Repeat_Count)) +
      geom_line(color = "black") +
      labs(x = "Position on Chromosome", y = "Number of Repeats", 
           title = paste("Number of Repeats for", chr)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # Print the plot to display it
    print(ggplot_object)
  }
}
