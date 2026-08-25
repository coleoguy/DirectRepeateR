# DirectRepeateR 1.0.1

## Bug fixes

* `GetRepeats()` no longer merges stale per-chromosome CSV files left over
  from previous runs. Intermediate files now go to a unique temporary
  directory by default (new `outdir` argument, threaded down to the C++
  routine); if a persistent `outdir` is supplied, pre-existing
  `*_condensed.csv` files in it are deleted before the run.
* The C++ FASTA reader now strips Windows-style carriage returns (CRLF
  line endings previously shifted every coordinate and corrupted
  chromosome names), uppercases sequence on read-in (soft-masked,
  lowercase genomes previously missed repeats silently), and skips query
  chunks containing `N` (assembly-gap N-runs previously matched each
  other and produced large false repeat blocks).
* Chromosome names are sanitized before being used as file names, so
  NCBI-style headers (e.g. `>gi|...|ref|NC_003279.8|`) no longer produce
  invalid paths; a failure to open an output file is now a hard error
  instead of a warning that silently dropped the chromosome.
* Chromosomes longer than 2^31 - 1 bases now raise an informative error
  instead of overflowing `int`; the union-find used to condense repeats
  is now iterative, so very large tandem arrays can no longer overflow
  the C stack.
* `PlotRepeats()` uses `linewidth` instead of the deprecated ggplot2
  `size` aesthetic and now (invisibly) returns the list of ggplot
  objects so plots can be saved or modified.
* `NAMESPACE` no longer imports knitr/rmarkdown (which are Suggests-only
  and made the package fail to install without them), no longer uses
  `exportPattern()`, and no longer exports the internal
  `run_combined_cpp()`.

## Changed behavior

* `ConvertToGFF()` now emits valid GFF3: nine columns (`seqid`, `source`,
  `type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`), the
  Sequence Ontology term `direct_repeat` for repeat copies (previously
  the invalid `repeat_copy`), and — when written to file — a
  `##gff-version 3` pragma with no column header row. The output file
  location is now controlled by a new `path` argument (default `NULL`,
  meaning no file is written) instead of the hard-coded
  `results/gff_output.txt`. The returned object is a nine-column data
  frame in GFF3 layout instead of the previous six-column table. The
  conversion is also vectorized (measured ~470x faster on the vignette
  data).

## Performance

* The core C++ scanner uses `std::string_view::find()` instead of a
  byte-at-a-time comparison loop (identical results, substantially
  faster on real chromosomes).
* `PlotRepeats()` bins repeats with `findInterval()`/`tabulate()` instead
  of an O(windows x repeats) `sapply()` loop (measured >100x faster on
  large inputs, identical counts).
* The FASTA reader accumulates sequence into a `std::string` directly,
  roughly halving transient memory on large genomes.

## Packaging

* `License` field corrected to `MIT + file LICENSE` with the standard
  two-line LICENSE template.
* The example C. elegans chromosome moved from `data/` to
  `inst/extdata/`; access it with
  `system.file("extdata", "C_elegans_chr1.fasta", package = "DirectRepeateR")`.
* Removed `LazyData` (the package has no `data/` datasets).
* Documentation corrected: `minlength` default is 50; `PlotRepeats()`
  argument is `data`; documented that tandem repeats with period shorter
  than `query_length` are reported as overlapping fragments.
* Added a testthat test suite.
