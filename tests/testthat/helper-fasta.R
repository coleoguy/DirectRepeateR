# Helpers for building small synthetic FASTA files with planted repeats.

# A 30-bp repeat unit whose 10-mers are unique (no internal periodicity).
repeat_unit <- "ACGGTTAACCGGATCCGGTACTGACTGGCA"

# Non-repetitive spacer (50 bp) and tail (20 bp): share no 10-mer with
# the repeat unit or each other.
spacer_50 <- "TCTCAGAGTGTGCACAGTCACGCGTATATGCGCTCTGTGACACTGAGAGA"
tail_20   <- "GATTACAGATTACAGGGCCC"

# chr1: repeat unit at positions 1-30 and again at 81-110 (130 bp total).
# With query_length = 10 the expected condensed row is
#   Start = 1, End = 30, Match_Start = 81, Match_End = 110.
chr1_seq <- paste0(repeat_unit, spacer_50, repeat_unit, tail_20)

# chr2: no repeats at all (70 bp of unique sequence).
chr2_seq <- paste0(spacer_50, tail_20)

# Write a FASTA file. `eol` lets tests exercise CRLF line endings.
write_fasta <- function(headers, seqs, path, width = 60, eol = "\n") {
  stopifnot(length(headers) == length(seqs))
  out <- character(0)
  for (i in seq_along(headers)) {
    out <- c(out, paste0(">", headers[i]))
    s <- seqs[i]
    starts <- seq(1, nchar(s), by = width)
    out <- c(out, substring(s, starts, pmin(starts + width - 1, nchar(s))))
  }
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(out, con, sep = eol)
  invisible(path)
}

# Run GetRepeats quietly (the C++ routine prints progress messages).
quiet_GetRepeats <- function(...) {
  res <- NULL
  capture.output(
    suppressMessages(res <- GetRepeats(...))
  )
  res
}
