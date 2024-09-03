ConvertToGFF <- function(data) {
  
  # Create a folder for results files if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Preallocate vectors to store GFF entries
  chromosomes <- vector("character", 3 * nrow(data))
  sources <- rep("GetRepeats", 3 * nrow(data))
  features <- vector("character", 3 * nrow(data))
  starts <- vector("integer", 3 * nrow(data))
  ends <- vector("integer", 3 * nrow(data))
  attributes <- vector("character", 3 * nrow(data))
  
  # Iterate over each row in data
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    
    # Generate ID for the repeat
    repeat_id <- paste0("Repeat_", i)
    
    # Index calculation
    idx <- (3 * i - 2):(3 * i)
    
    # Full length of the repeat (Start_Position of first copy to End_Position of second copy)
    chromosomes[idx[1]] <- row$Chromosome
    features[idx[1]] <- "repeat_region"
    starts[idx[1]] <- row$Start_Position
    ends[idx[1]] <- row$Match_End_Position
    attributes[idx[1]] <- paste0("ID=", repeat_id)
    
    # First copy of the repeat
    chromosomes[idx[2]] <- row$Chromosome
    features[idx[2]] <- "repeat_copy"
    starts[idx[2]] <- row$Start_Position
    ends[idx[2]] <- row$End_Position
    attributes[idx[2]] <- paste0("Parent=", repeat_id, ";Copy=1")
    
    # Second copy of the repeat
    chromosomes[idx[3]] <- row$Chromosome
    features[idx[3]] <- "repeat_copy"
    starts[idx[3]] <- row$Match_Position
    ends[idx[3]] <- row$Match_End_Position
    attributes[idx[3]] <- paste0("Parent=", repeat_id, ";Copy=2")
  }
  
  # Combine vectors into a data frame
  gff_output <- data.frame(Chromosome = chromosomes,
                           Source = sources,
                           Feature = features,
                           Start = starts,
                           End = ends,
                           Attributes = attributes,
                           stringsAsFactors = FALSE)
  
  # Write the GFF output to a file
  write.table(gff_output, file = "results/gff_output.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

