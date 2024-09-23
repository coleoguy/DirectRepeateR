GetRepeats <- function(file, query_length = NULL, maxdist = NULL, minlength = NULL) {
  
  # Set default values if parameters are NULL
  if (is.null(query_length)) {
    query_length <- 25  # Default length for query sequences
  }
  if (is.null(maxdist)) {
    maxdist <- 20000  # Default maximum distance between repeat copies
  }
  if (is.null(minlength)) {
    minlength <- 25  # Default minimum length of repeats to consider
  }
  
  # Read in FASTA file, obtain chromosome lengths and names
  genome <- read.fasta(file)  # Load the genome sequences from the FASTA file
  SeqLength <- sapply(genome, length)  # Determine the length of each chromosome
  chromosome_names <- names(genome)  # Get the names of the chromosomes
  
  # Function to create an empty results data.table with a specified number of rows
  create_empty_results <- function(n = 1000000) {
    data.frame(Chromosome = character(n),  # Placeholder for chromosome names
               Start_Position = integer(n),  # Placeholder for start positions of repeats
               Match_Position = integer(n),  # Placeholder for start positions of repeat matches
               End_Position = integer(n),  # Placeholder for end positions of repeats
               Match_End_Position = integer(n))  # Placeholder for end positions of repeat matches
  }
  
  # Create a folder for temporary files if it doesn't exist
  if (!dir.exists("temp")) {
    dir.create("temp")
  }
  
  # Create a folder for results files if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Function to condense contiguous results into single entries
  condense_results <- function(data) {
    chr <- unique(data$Chromosome)  # Get unique chromosome names
    condensed_results <- list()  # Initialize a list to store condensed results
    
    # Iterate across results for each chromosome
    for (i in seq_along(chr)) {
      dat_chr <- data[data$Chromosome == chr[i],]  # Subset data to the current chromosome
      results_chr <- list()  # Initialize a list to store results for the current chromosome
      result_index <- 1  # Initialize an index for storing results
      
      j <- 1 
      n <- nrow(dat_chr)  # Get the number of rows for the current chromosome
      
      # Loop through each row in the chromosome data
      while (j <= n) {
        start_pos <- dat_chr$Start_Position[j]  # Start position of the repeat
        match_start_pos <- dat_chr$Match_Position[j]  # Start position of the repeat match
        end_pos <- dat_chr$End_Position[j]  # End position of the repeat
        match_end_pos <- dat_chr$Match_End_Position[j]  # End position of the repeat match
        
        curset <- j  # Initialize current set with the current row
        
        k <- j + 1  # Move to the next row
        # While loop to check for contiguous repeat sequences
        while (k <= n && dat_chr$Start_Position[k] <= (end_pos + 1)) {
          # If the start position and match position of the next row are contiguous
          if (dat_chr$Start_Position[k] == (end_pos + 1) && 
              dat_chr$Match_Position[k] == (match_end_pos + 1)) {
            end_pos <- dat_chr$End_Position[k]  # Update end position
            match_end_pos <- dat_chr$Match_End_Position[k]  # Update match end position
            curset <- c(curset, k)  # Add the current row to the set
          }
          k <- k + 1  # Move to the next row
        }
        
        # Create a condensed result entry for the current set
        condensed_set <- data.table(Chromosome = chr[i], 
                                    Start_Position = start_pos, 
                                    End_Position = end_pos, 
                                    Match_Position = match_start_pos, 
                                    Match_End_Position = match_end_pos)
        results_chr[[result_index]] <- condensed_set  # Store the condensed result
        result_index <- result_index + 1  # Increment the result index
        
        j <- curset[length(curset)] + 1  # Move to the next unprocessed row
      }
      
      condensed_results[[chr[i]]] <- rbindlist(results_chr)  # Combine all results for the chromosome
    }
    
    return(rbindlist(condensed_results))  # Return the combined results for all chromosomes
  }
  
  # Initialize variables
  raw_results <- create_empty_results()  # Create an empty data frame to store raw results
  insert_idx <- 1  # Initialize index for inserting results into the data frame
  file_count <- 1  # Initialize a counter for temporary files
  
  # Loop through each chromosome in the genome
  for (i in seq_along(genome)) {
    sequence <- paste(genome[[i]], collapse = "")  # Convert the chromosome sequence to a single string
    seq_length <- SeqLength[i]  # Get the length of the chromosome
    starts <- seq(from = 1, by = query_length, length.out = floor(seq_length / query_length))  # Generate start positions based on query length
    
    # Loop through each start position to search for repeats
    for (j in seq_along(starts)) {
      start_pos <- starts[j]  # Get the current start position
      upper_bound <- min(start_pos + query_length - 1 + maxdist, seq_length)  # Define the upper bound for searching
      pattern <- substr(sequence, start_pos, start_pos + query_length - 1)  # Extract the query pattern from the sequence
      search_field <- substr(sequence, start_pos + query_length, upper_bound)  # Define the search field
      
      # Search for matches of the query pattern in the search field
      start_matches <- gregexpr(pattern = pattern, text = search_field, fixed = TRUE)[[1]]
      
      # If matches are found
      if (any(start_matches > 0)) {
        curr_match_pos <- start_pos + query_length - 1 + start_matches  # Calculate the match positions
        num_matches <- length(curr_match_pos)  # Get the number of matches
        
        # Check if the results data frame is full
        if ((insert_idx + num_matches - 1) > nrow(raw_results)) {
          # If full, remove empty rows and condense the results
          temp_results <- raw_results[1:(insert_idx - 1), ]
          temp_results <- temp_results[temp_results$Chromosome != "", ]
          temp_results <- condense_results(temp_results)
          
          # Write the condensed results to a temporary file
          fwrite(temp_results, paste0("temp/temp.", file_count, ".csv"))
          file_count <- file_count + 1  # Increment the file counter
          raw_results <- create_empty_results()  # Create a new empty results data frame
          insert_idx <- 1  # Reset the insert index
        }
        
        # Add the found matches to the results data frame
        raw_results[insert_idx:(insert_idx + num_matches - 1), ] <- data.frame(
          Chromosome = rep(chromosome_names[i], num_matches),  # Repeat the chromosome name for each match
          Start_Position = rep(start_pos, num_matches),  # Repeat the start position for each match
          Match_Position = curr_match_pos,  # Match positions
          End_Position = rep(start_pos + query_length - 1, num_matches),  # End positions
          Match_End_Position = curr_match_pos + query_length - 1  # Match end positions
        )
        
        # Update the insert index
        insert_idx <- insert_idx + num_matches
      }
    }
  }
  
  # Write any remaining results to a file
  if (insert_idx > 1) {
    temp_results <- raw_results[1:(insert_idx - 1), ]  # Extract the filled portion of the results data frame
    temp_results <- temp_results[temp_results$Chromosome != "", ]  # Remove rows with empty chromosome names
    temp_results <- condense_results(temp_results)  # Condense the results
    fwrite(temp_results, paste0("temp/temp.", file_count, ".csv"))  # Write to a temporary file
  }
  
  # Combine all temporary files into one data frame
  temp_files <- list.files("temp", pattern = "temp.*\\.csv", full.names = TRUE)  # List all temp files
  combined_results <- rbindlist(lapply(temp_files, fread))  # Read and combine all temp files
  
  # Remove the temporary files
  file.remove(temp_files)
  
  # Apply the condensing function to the combined results
  condensed_results <- condense_results(combined_results)
  
  # Calculate the length of each repeat
  condensed_results[, Repeat_Length := End_Position - Start_Position]
  
  # Filter out repeats shorter than the minimum length
  filtered_results <- condensed_results[Repeat_Length > minlength]
  
  # Remove the Repeat_Length column before saving
  final_output <- filtered_results[, .(Chromosome, Start_Position, End_Position, Match_Position, Match_End_Position)]
  
  # Save the final results to a CSV file
  fwrite(final_output, "results/condensed_results.csv")
  
  return(final_output)  # Return the final results
}
