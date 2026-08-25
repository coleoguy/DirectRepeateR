\name{GetRepeats}
\alias{GetRepeats}
\title{Detection of Direct Repeats in Genome Assemblies}
\description{
  The \code{GetRepeats} function calls a compiled C++ routine (via \pkg{Rcpp}) 
  to detect direct repeats in a FASTA file. It stores chromosomal CSV  
  files in an intermediate directory (a unique temporary directory by 
  default), then merges and filters them, returning a combined table of 
  repeat locations.
}
\usage{
GetRepeats(
  file,
  query_length = 25,
  maxdist      = 20000,
  minlength    = 50,
  merged_csv   = "merged_results.csv",
  outdir       = tempfile("chromosome_results_")
)
}
\arguments{
  \item{file}{
    A character string specifying the path to the genome assembly file 
    in FASTA format.
  }
  \item{query_length}{
    An integer specifying the length of the query chunks used by the 
    direct-repeat search algorithm. Default is \code{25}.
  }
  \item{maxdist}{
    An integer specifying how far beyond each chunk the function will 
    search for repeats. Default is \code{20000}.
  }
  \item{minlength}{
    An integer specifying the minimum length of repeats to keep.
    Repeats shorter than this are discarded. Default is \code{50}.
  }
  \item{merged_csv}{
    A character string specifying the path of a final merged CSV file. 
    Defaults to \code{"merged_results.csv"} in the working directory. 
    If \code{NULL}, no merged CSV is saved to disk.
  }
  \item{outdir}{
    A character string specifying the directory for the intermediate 
    per-chromosome CSV files. Defaults to a unique temporary directory 
    so results can never be contaminated by files left over from 
    previous runs. If a persistent directory is supplied, any 
    pre-existing \code{"_condensed.csv"} files in it are deleted before 
    the run.
  }
}
\details{
  The search uppercases sequence on read-in (so soft-masked, lowercase 
  genomes are handled), strips Windows-style carriage returns, and skips 
  query chunks containing \code{N} (assembly gaps).

  Tandem (periodic) repeats whose period is shorter than 
  \code{query_length} are reported as multiple overlapping fragments 
  rather than a single merged interval; overlapping rows in the output 
  can be merged post hoc if a single interval per tandem array is 
  desired.
}
\value{
  A \code{data.table} with the following columns:
  \describe{
    \item{Chromosome}{Chromosome/sequence the repeat is on.}
    \item{Start}{Start position of the first occurrence of the repeat.}
    \item{End}{End position of the first occurrence of the repeat.}
    \item{Match_Start}{Start position of the matching repeat copy.}
    \item{Match_End}{End position of the matching repeat copy.}
  }
}
\author{
  Megan Copeland, Andres Barboza, Heath Blackmon
}
\references{
  \url{https://github.com/coleoguy/DirectRepeateR}
}
\examples{
\dontrun{
  # Example usage:
  final_dt <- GetRepeats(
    file         = "data/human_filt.fa",
    query_length = 25,
    maxdist      = 20000,
    minlength    = 50,
    merged_csv   = "my_final_merged.csv"
  )
}
}
