# Tests for ConvertToGFF structure and GFF3 validity.

toy_repeats <- data.frame(
  Chromosome  = c("chr1", "chr1", "chr2"),
  Start       = c(100L, 5000L, 300L),
  End         = c(150L, 5100L, 400L),
  Match_Start = c(500L, 9000L, 700L),
  Match_End   = c(550L, 9100L, 750L),
  stringsAsFactors = FALSE
)

test_that("ConvertToGFF returns nine GFF3 columns, three rows per repeat", {
  gff <- ConvertToGFF(toy_repeats)
  expect_equal(names(gff),
               c("seqid", "source", "type", "start", "end",
                 "score", "strand", "phase", "attributes"))
  expect_equal(nrow(gff), 3L * nrow(toy_repeats))
  expect_true(all(gff$type %in% c("repeat_region", "direct_repeat")))
  expect_true(all(gff$score == "."))
  expect_true(all(gff$strand == "."))
  expect_true(all(gff$phase == "."))
})

test_that("ConvertToGFF coordinates are correct for every repeat row", {
  gff <- ConvertToGFF(toy_repeats)
  for (i in seq_len(nrow(toy_repeats))) {
    block <- gff[(3 * (i - 1) + 1):(3 * i), ]
    r <- toy_repeats[i, ]
    # repeat_region spans first copy start to second copy end
    expect_equal(block$type[1], "repeat_region")
    expect_equal(block$start[1], r$Start)
    expect_equal(block$end[1],   r$Match_End)
    expect_equal(block$seqid,    rep(r$Chromosome, 3))
    # copy 1
    expect_equal(block$start[2], r$Start)
    expect_equal(block$end[2],   r$End)
    # copy 2
    expect_equal(block$start[3], r$Match_Start)
    expect_equal(block$end[3],   r$Match_End)
    # Parent linkage
    id <- sub("^ID=([^;]+).*$", "\\1", block$attributes[1])
    expect_match(block$attributes[2], paste0("Parent=", id), fixed = TRUE)
    expect_match(block$attributes[3], paste0("Parent=", id), fixed = TRUE)
  }
})

test_that("ConvertToGFF writes a valid GFF3 file when path is given", {
  out <- tempfile(fileext = ".gff3")
  gff <- ConvertToGFF(toy_repeats, path = out)
  expect_true(file.exists(out))
  lines <- readLines(out)
  expect_equal(lines[1], "##gff-version 3")
  # No header row; one line per feature plus the pragma
  expect_equal(length(lines), 1L + nrow(gff))
  # Every feature line has exactly 9 tab-separated fields
  fields <- strsplit(lines[-1], "\t", fixed = TRUE)
  expect_true(all(lengths(fields) == 9L))
})

test_that("ConvertToGFF writes no file by default and handles empty input", {
  wd_before <- list.files()
  gff <- ConvertToGFF(toy_repeats)
  expect_equal(list.files(), wd_before)  # no results/ dir created
  
  empty <- toy_repeats[0, ]
  gff0 <- ConvertToGFF(empty)
  expect_equal(nrow(gff0), 0L)
  expect_equal(names(gff0),
               c("seqid", "source", "type", "start", "end",
                 "score", "strand", "phase", "attributes"))
})
