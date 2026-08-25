# Regression tests for GetRepeats and the underlying C++ scanner.

expected_chr1 <- function(dt) {
  # One condensed repeat: positions 1-30 matching 81-110
  dt <- dt[dt$Chromosome == "chr1", ]
  expect_equal(nrow(dt), 1L)
  expect_equal(dt$Start, 1L)
  expect_equal(dt$End, 30L)
  expect_equal(dt$Match_Start, 81L)
  expect_equal(dt$Match_End, 110L)
}

test_that("GetRepeats finds a planted repeat in a multi-chromosome FASTA", {
  fa <- tempfile(fileext = ".fasta")
  write_fasta(c("chr1", "chr2"), c(chr1_seq, chr2_seq), fa)
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 30, merged_csv = NULL)
  expected_chr1(res)
  # chr2 has no repeats and must not appear
  expect_false("chr2" %in% res$Chromosome)
})

test_that("CRLF (Windows) line endings give identical results (drpt-02)", {
  fa <- tempfile(fileext = ".fasta")
  write_fasta(c("chr1", "chr2"), c(chr1_seq, chr2_seq), fa, eol = "\r\n")
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 30, merged_csv = NULL)
  expected_chr1(res)
  # Chromosome name must not carry a trailing \r
  expect_true(all(res$Chromosome == "chr1"))
})

test_that("lowercase soft-masked sequence is matched (drpt-03)", {
  fa <- tempfile(fileext = ".fasta")
  # Soft-mask the second copy of the repeat only
  chr1_soft <- paste0(repeat_unit, spacer_50, tolower(repeat_unit), tail_20)
  write_fasta("chr1", chr1_soft, fa)
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 30, merged_csv = NULL)
  expected_chr1(res)
})

test_that("N-runs do not generate false repeats (drpt-03)", {
  fa <- tempfile(fileext = ".fasta")
  # Two 40-bp N-runs separated by unique sequence: naive matching would
  # report them as a direct repeat.
  n_run <- strrep("N", 40)
  chrN <- paste0(n_run, spacer_50, n_run, tail_20)
  write_fasta("chrN", chrN, fa)
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 10, merged_csv = NULL)
  expect_equal(nrow(res), 0L)
})

test_that("stale condensed CSVs from previous runs are not merged (drpt-01)", {
  outdir <- tempfile("chromosome_results_")
  dir.create(outdir)
  # Plant a stale file as if from a previous run on another genome
  stale <- data.frame(Start = 999L, End = 1200L,
                      Match_Start = 5000L, Match_End = 5201L)
  write.csv(stale, file.path(outdir, "stale_chrom_condensed.csv"),
            row.names = FALSE)
  
  fa <- tempfile(fileext = ".fasta")
  write_fasta("chr1", chr1_seq, fa)
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 30, merged_csv = NULL,
                          outdir = outdir)
  expect_false("stale_chrom" %in% res$Chromosome)
  expected_chr1(res)
})

test_that("special characters in chromosome names are sanitized (drpt-07)", {
  fa <- tempfile(fileext = ".fasta")
  write_fasta("gi|123|ref|NC_000001.1|", chr1_seq, fa)
  res <- quiet_GetRepeats(fa, query_length = 10, maxdist = 200,
                          minlength = 30, merged_csv = NULL)
  # The chromosome must not be silently dropped
  expect_equal(nrow(res), 1L)
  expect_equal(res$Chromosome, "gi_123_ref_NC_000001.1_")
})
