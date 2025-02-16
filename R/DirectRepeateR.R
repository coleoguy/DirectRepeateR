#' Get Direct Repeats
#'
#' This function:
#'   1. Calls a C++ routine to search for direct repeats in a FASTA file,
#'      always storing per-chromosome CSVs in "chromosome_results".
#'   2. Merges those "_condensed.csv" files into one table.
#'   3. Filters out repeats below a specified minimum length.
#'   4. By default, saves the merged results to "merged_results.csv"
#'      in the working directory, unless you specify `merged_csv = NULL`.
#'
#' @param file         Path to the input FASTA file.
#' @param query_length Size of the chunks used for searching. Default is 25.
#' @param maxdist      Maximum distance beyond each chunk to search. Default is 20000.
#' @param minlength    Minimum repeat length to retain in the final table. Default is 25.
#' @param merged_csv   A path/filename for the merged CSV. Defaults to 
#'   `"merged_results.csv"` in the working directory. If \code{NULL}, 
#'   no merged file is written.
#'
#' @return A `data.table` with columns:
#'   \itemize{
#'     \item Start_Position
#'     \item End_Position
#'     \item Match_Position
#'     \item Match_End_Position
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # By default, this saves "merged_results.csv" to your working directory:
#'   final_dt <- GetRepeats("data/human_filt.fa")
#'
#'   # If you want NO file output, set merged_csv = NULL:
#'   final_dt2 <- GetRepeats("data/human_filt.fa", merged_csv = NULL)
#'
#'   # Or provide a custom path:
#'   final_dt3 <- GetRepeats("data/human_filt.fa", merged_csv = "my_output.csv")
#' }
GetRepeats <- function(query_length = 25,
                       maxdist      = 20000,
                       minlength    = 50,
                       merged_csv   = "merged_results.csv") {
  
  library(data.table)  # Use data.table for efficiency
  
  # 1) Run the C++ routine to generate per-chromosome CSVs
  run_combined_cpp(
    fasta_path   = "your_fasta_file.fa",  # Replace with your actual FASTA file path
    query_length = query_length,
    maxdist      = maxdist
  )

  # 2) Collect all "_condensed.csv" files from "chromosome_results"
  outdir <- "chromosome_results"
  condensed_files <- list.files(
    path       = outdir,
    pattern    = "_condensed\\.csv$",
    full.names = TRUE
  )

  # Debug: Check if files are found
  if (length(condensed_files) == 0) {
    stop("No condensed CSV files found in '", outdir, "'. Possibly no repeats identified.")
  }
  
  # 3) Function to process each file
  process_csv <- function(csv_file) {
    # Read CSV
    df <- read.csv(csv_file, stringsAsFactors = FALSE)
    
    # Skip empty files
    if (nrow(df) == 0) {
      print(paste("Skipping empty file:", csv_file))
      return(NULL)
    }
    
    # Extract chromosome name from filename
    chrom_name <- sub("_condensed\\.csv$", "", basename(csv_file))
    
    # Add Chromosome column as the first column
    df$Chromosome <- chrom_name
    df <- df[, c(5, 1:4)]  # Reorder columns

    # Sort by column 2 (Start_Position) then column 4 (Match_Position)
    df <- df[order(df[[2]], df[[4]]), ]
    
    # Compute repeat length
    df$Repeat_Length <- (df$End_Position - df$Start_Position) + 1
    
    # Apply filtering (keeping only repeats >= minlength)
    df <- df[df$Repeat_Length >= minlength, ]
    
    # Remove temporary Repeat_Length column
    df$Repeat_Length <- NULL

    return(df)
  }

  # 4) Process all files
  list_of_dfs <- lapply(condensed_files, process_csv)

  # Remove NULL elements (empty files)
  list_of_dfs <- Filter(Negate(is.null), list_of_dfs)

  # 5) Merge all filtered data frames
  merged_df <- do.call(rbind, list_of_dfs)

  # Debug: Print merged table preview
  print("Merged data preview:")
  print(head(merged_df))

  # 6) Save merged results
  write.csv(merged_df, merged_csv, row.names = FALSE)

  # Debug: Confirm file saved
  print(paste("Merged file saved as:", merged_csv))

  return(merged_df)
}

# Run the function
merged_results <- GetRepeats()
