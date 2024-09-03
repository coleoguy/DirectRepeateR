\name{GetRepeats}
\alias{GetRepeats}
\title{Detection of Direct Repeats in Genome Assemblies}
\description{
  The \code{GetRepeats} function searches for and stores the locations 
  of direct repeats within genome assemblies provided in FASTA format.
}
\usage{
GetRepeats(file, query_length = NULL, maxdist = NULL, minlength = NULL)
}
\arguments{
  \item{file}{
    A character string specifying the path to the genome assembly file 
    in FASTA format. The function uses the \pkg{seqinr} package to read 
    the FASTA file.
  }
  \item{query_length}{
    An integer specifying the length of the query sequence used to search 
    for direct repeats. The genome sequence is segmented into subsequences 
    of this length to identify repeated patterns. Default is \code{25}.
  }
  \item{maxdist}{
    An integer specifying the maximum distance between repeat copies 
    within which the function will search for direct repeats. This sets 
    a limit on how far apart repeat sequences can be to be considered 
    related. Default is \code{20000}.
  }
  \item{minlength}{
    An integer specifying the minimum length of sequences to be considered 
    as direct repeats. Only repeats of this length or longer will be detected.
    Default is \code{25}.
  }
}
\value{
  Returns a data frame with five columns detailing:
  \describe{
    \item{Chromosome}{Chromosome identifier where the repeat is located.}
    \item{Start_Position}{Start position of the first occurrence of the repeat.}
    \item{End_Position}{End position of the first occurrence of the repeat.}
    \item{Match_Position}{Start position of the second occurrence of the repeat 
    (the repeat copy).}
    \item{Match_End_Position}{End position of the second occurrence of the repeat 
    (the repeat copy).}
  }
}
\author{
  Megan Copeland, Heath Blackmon
}
\references{
  \url{https://github.com/coleoguy/DirectRepeateR}
}
\examples{
# The following code is an example, but it is commented out to prevent execution.
# GetRepeats("genome.fa", query_length = 25, maxdist = 20000, minlength = 50)
}
\keyword{direct repeats}
