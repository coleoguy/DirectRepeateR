# Tests for PlotRepeats binning counts and return value.

# Reference (old) implementation of the binning
ref_counts <- function(mid, windows, window_size) {
  sapply(windows, function(w) {
    sum(mid >= w & mid < (w + window_size))
  })
}

make_repeat_df <- function(mids, chrom = "chr1") {
  # Build a repeat table whose midpoints equal `mids` exactly:
  # place both copies as zero-length-ish intervals centered on mid.
  data.frame(
    Chromosome  = chrom,
    Start       = mids - 2,
    End         = mids,
    Match_Start = mids,
    Match_End   = mids + 2,
    stringsAsFactors = FALSE
  )
}

test_that("PlotRepeats bin counts match the reference implementation", {
  set.seed(42)
  mids <- sort(sample(10:99990, 500))
  df <- make_repeat_df(mids)
  window_size <- 10000
  step_size <- 10000
  
  chr_length <- max(c(df$End, df$Match_End))
  windows <- seq(1, chr_length, by = step_size)
  mid <- ((df$Start + df$End) / 2 + (df$Match_Start + df$Match_End) / 2) / 2
  expected <- ref_counts(mid, windows, window_size)
  
  pdf(NULL)  # swallow plot output
  on.exit(dev.off())
  plots <- PlotRepeats(df, window_size = window_size, step_size = step_size)
  
  built <- ggplot2::ggplot_build(plots[["chr1"]])
  counts <- built$plot$data$Repeat_Count
  expect_equal(counts, expected)
  expect_equal(sum(counts), length(mids))
})

test_that("PlotRepeats handles overlapping windows (step < window)", {
  set.seed(7)
  mids <- sort(sample(10:49990, 200))
  df <- make_repeat_df(mids)
  window_size <- 10000
  step_size <- 5000
  
  chr_length <- max(c(df$End, df$Match_End))
  windows <- seq(1, chr_length, by = step_size)
  mid <- ((df$Start + df$End) / 2 + (df$Match_Start + df$Match_End) / 2) / 2
  expected <- ref_counts(mid, windows, window_size)
  
  pdf(NULL)
  on.exit(dev.off())
  plots <- PlotRepeats(df, window_size = window_size, step_size = step_size)
  counts <- plots[["chr1"]]$data$Repeat_Count
  expect_equal(counts, expected)
})

test_that("PlotRepeats returns one ggplot per chromosome, invisibly", {
  df <- rbind(make_repeat_df(c(1000, 2000, 3000), "chr1"),
              make_repeat_df(c(1500, 2500), "chr2"))
  pdf(NULL)
  on.exit(dev.off())
  res <- withVisible(PlotRepeats(df))
  expect_false(res$visible)
  expect_named(res$value, c("chr1", "chr2"))
  expect_s3_class(res$value[["chr1"]], "ggplot")
  expect_s3_class(res$value[["chr2"]], "ggplot")
})
