ConvertToGFF <- function(data) {
  
  # Create a folder for results files if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Initialize an empty list to store GFF entries
  gff_list <- list()
  
  # Iterate over each row in the data
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    
    # Generate ID for the repeat
    repeat_id <- paste0("Repeat_", i)
    
    # Full length of the repeat (Start_Position of first copy to Match_End_Position of second copy)
    gff_list[[length(gff_list) + 1]] <- data.frame(
      Chromosome = row$Chromosome,
      Source = "GetRepeats",
      Feature = "repeat_region",
      Start = row$Start_Position,
      End = row$Match_End_Position,
      Attributes = paste0("ID=", repeat_id),
      stringsAsFactors = FALSE
    )
    
    # First copy of the repeat
    gff_list[[length(gff_list) + 1]] <- data.frame(
      Chromosome = row$Chromosome,
      Source = "GetRepeats",
      Feature = "repeat_copy",
      Start = row$Start_Position,
      End = row$End_Position,
      Attributes = paste0("Parent=", repeat_id, ";Copy=1"),
      stringsAsFactors = FALSE
    )
    
    # Second copy of the repeat
    gff_list[[length(gff_list) + 1]] <- data.frame(
      Chromosome = row$Chromosome,
      Source = "GetRepeats",
      Feature = "repeat_copy",
      Start = row$Match_Position,
      End = row$Match_End_Position,
      Attributes = paste0("Parent=", repeat_id, ";Copy=2"),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine list into a data frame
  gff_output <- do.call(rbind, gff_list)
  
  # Write the GFF output to a file
  write.table(gff_output, file = "results/gff_output.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
  
  return(gff_output)
}
