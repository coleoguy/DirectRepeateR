ConvertToGFF <- function(data, path = NULL) {
  
  n <- nrow(data)
  
  if (n == 0) {
    gff_output <- data.frame(
      seqid      = character(0),
      source     = character(0),
      type       = character(0),
      start      = integer(0),
      end        = integer(0),
      score      = character(0),
      strand     = character(0),
      phase      = character(0),
      attributes = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    repeat_id <- paste0("Repeat_", seq_len(n))
    
    # Three GFF rows per repeat, interleaved in the order:
    # repeat_region (full span), then the two direct_repeat copies.
    idx <- rep(seq_len(n), each = 3L)
    which_row <- rep.int(1:3, n)
    
    starts <- ifelse(which_row == 3L, data$Match_Start[idx], data$Start[idx])
    ends   <- ifelse(which_row == 2L, data$End[idx], data$Match_End[idx])
    
    types <- c("repeat_region", "direct_repeat", "direct_repeat")[which_row]
    
    attributes <- character(3L * n)
    attributes[which_row == 1L] <- paste0("ID=", repeat_id)
    attributes[which_row == 2L] <- paste0("ID=", repeat_id, ".copy1;Parent=",
                                          repeat_id)
    attributes[which_row == 3L] <- paste0("ID=", repeat_id, ".copy2;Parent=",
                                          repeat_id)
    
    gff_output <- data.frame(
      seqid      = as.character(data$Chromosome)[idx],
      source     = "DirectRepeateR",
      type       = types,
      start      = starts,
      end        = ends,
      score      = ".",
      strand     = ".",
      phase      = ".",
      attributes = attributes,
      stringsAsFactors = FALSE
    )
  }
  
  # Optionally write a valid GFF3 file (pragma line, no column header)
  if (!is.null(path)) {
    con <- file(path, open = "wt")
    on.exit(close(con))
    writeLines("##gff-version 3", con)
    write.table(gff_output, file = con, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
  }
  
  return(gff_output)
}
