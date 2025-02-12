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
GetRepeats <- function(file,
                       query_length = 25,
                       maxdist      = 20000,
                       minlength    = 25,
                       merged_csv   = "merged_results.csv") {
  
  # 1) Calls the C++ routine, which writes per-chromosome CSV files
  #    into the folder "chromosome_results". This folder is not configurable.
  run_combined_cpp(
    fasta_path   = file,
    query_length = query_length,
    maxdist      = maxdist
  )
  
  # 2) Collect the new "_condensed.csv" files from "chromosome_results"
  outdir <- "chromosome_results"
  condensed_files <- list.files(
    path       = outdir,
    pattern    = "_condensed\\.csv$",
    full.names = TRUE
  )
  
  if (length(condensed_files) == 0) {
    message("No condensed CSV files found in '", outdir, "'. Possibly no repeats identified.")
    return(data.table::data.table())
  }
  
  # 3) Read files, compute length for filtering, then remove it
  list_of_dt <- lapply(condensed_files, function(csv_file) {
    dt <- data.table::fread(csv_file)
    data.table::setorder(dt, Start_Position, Match_Position)
    dt[, temp_length := (End_Position - Start_Position) + 1]  # for filtering
    dt <- dt[temp_length >= minlength]
    dt[, temp_length := NULL]  # remove before returning
    dt
  })
  
  # 4) Merge into one data.table
  merged_dt <- data.table::rbindlist(list_of_dt, use.names = TRUE, fill = TRUE)
  
  # 5) By default, save "merged_results.csv" unless merged_csv is NULL
  if (!is.null(merged_csv)) {
    data.table::fwrite(merged_dt, merged_csv)
  }
  
  # 6) Return the final data.table
  return(merged_dt)
}
