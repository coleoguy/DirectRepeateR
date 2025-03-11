\name{GetRepeats}
\alias{GetRepeats}
\title{Detection of Direct Repeats in Genome Assemblies}
\description{
  The \code{GetRepeats} function calls a compiled C++ routine (via \pkg{Rcpp}) 
  to detect direct repeats in a FASTA file. It stores chromosomal CSV  
  files in the directory \code{"chromosome_results"}, then merges 
  and filters them, returning a combined table of repeat locations.
}
\usage{
GetRepeats(
  file,
  query_length = 25,
  maxdist      = 20000,
  minlength    = 25,
  merged_csv   = merged_results.csv
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
    Repeats shorter than this are discarded. Default is \code{25}.
  }
  \item{merged_csv}{
    (Optional) A character string specifying the path of a final merged 
    CSV file. If \code{NULL} (default), no merged CSV is saved to disk.
  }
}
\value{
  A \code{data.table} with the following columns:
  \describe{
    \item{Start}{Start position of the first occurrence of the repeat.}
    \item{End}{End position of the first occurrence of the repeat.}
    \item{Match_Start}{Start position of the matching repeat copy.}
    \item{Match_End}{End position of the matching repeat copy.}
  }
  
  The chromosomal intermediate CSV files are written to a folder named 
  \code{"chromosome_results"} in the working directory.
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
