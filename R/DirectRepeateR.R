#' Get Direct Repeats
#'
#' This function:
#'   1. Calls a C++ routine to search for direct repeats in a FASTA file,
#'      storing per-chromosome CSVs in an intermediate directory
#'      (a unique temporary directory by default).
#'   2. Merges those "_condensed.csv" files into one table.
#'   3. Filters out repeats below a specified minimum length.
#'   4. By default, saves the merged results to "merged_results.csv"
#'      in the working directory, unless you specify `merged_csv = NULL`.
#'
#' @param file         Path to the input FASTA file.
#' @param query_length Size of the chunks used for searching. Default is 25.
#' @param maxdist      Maximum distance beyond each chunk to search. Default is 20000.
#' @param minlength    Minimum repeat length to retain in the final table. Default is 50.
#' @param merged_csv   A path/filename for the merged CSV. Defaults to 
#'   `"merged_results.csv"` in the working directory. If \code{NULL}, 
#'   no merged file is written.
#' @param outdir       Directory for the intermediate per-chromosome CSV
#'   files. Defaults to a unique temporary directory so results can never
#'   be contaminated by files left over from previous runs. If you supply
#'   a persistent directory, any pre-existing "_condensed.csv" files in it
#'   are deleted before the run.
#'
#' @return A `data.table` with columns:
#'   \itemize{
#'     \item Chromosome
#'     \item Start
#'     \item End
#'     \item Match_Start
#'     \item Match_End
#'   }
#'
#' @details
#' Tandem (periodic) repeats whose period is shorter than
#' \code{query_length} are reported as multiple overlapping fragments
#' rather than one merged interval; overlapping rows in the output can be
#' merged post hoc if a single interval per array is desired.
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
GetRepeats <- function(file,
                       query_length = 25,
                       maxdist      = 20000,
                       minlength    = 50,
                       merged_csv   = "merged_results.csv",
                       outdir       = tempfile("chromosome_results_")) {
  
  # Defensively remove any stale per-chromosome files from a previous run
  # so they cannot be merged into the new results.
  if (dir.exists(outdir)) {
    stale <- list.files(outdir, pattern = "_condensed\\.csv$", full.names = TRUE)
    if (length(stale) > 0) {
      unlink(stale)
    }
  }
  
  # 1) Calls the C++ routine, which writes per-chromosome CSV files
  run_combined_cpp(
    fasta_path   = file,  # Use file argument instead of hardcoding
    query_length = query_length,
    maxdist      = maxdist,
    outdir       = outdir
  )
  
  # 2) Collect the new "_condensed.csv" files from the intermediate dir
  condensed_files <- list.files(
    path       = outdir,
    pattern    = "_condensed\\.csv$",
    full.names = TRUE
  )

  if (length(condensed_files) == 0) {
    message("No condensed CSV files found in '", outdir, "'. Possibly no repeats identified.")
    return(data.table::data.table())
  }
  
  # 3) Read files, extract chromosome name, sort, and filter
  list_of_dt <- lapply(condensed_files, function(csv_file) {
    dt <- data.table::fread(csv_file)
    
    # Extract chromosome name
    chrom_name <- sub("_condensed\\.csv$", "", basename(csv_file)) 

    # Add Chromosome column as the first column
    dt[, Chromosome := chrom_name]
    setcolorder(dt, c("Chromosome", setdiff(names(dt), "Chromosome")))

    # Sort by column 2 (Start) then column 4 (Match_Start)
    setorder(dt, Start, Match_Start)

    # Compute length for filtering
    dt[, temp_length := (End - Start) + 1]

    # Filter out short repeats
    dt <- dt[temp_length >= minlength]

    # Remove temporary column
    dt[, temp_length := NULL]

    return(dt)
  })

  # 4) Merge into one data.table
  merged_dt <- data.table::rbindlist(list_of_dt, use.names = TRUE, fill = TRUE)

  # 5) Save merged results if requested
  if (!is.null(merged_csv)) {
    data.table::fwrite(merged_dt, merged_csv)
  }

  # 6) Return the final data.table
  return(merged_dt)
}
