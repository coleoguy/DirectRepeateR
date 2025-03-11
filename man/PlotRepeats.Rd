\name{PlotRepeats}
\alias{PlotRepeats}
\title{Visualization of Rolling-Average Repeat Counts Across Chromosomes}
\description{
  The \code{PlotRepeats} function visualizes the distribution of direct 
  repeats across chromosomes by generating line plots 
  showing the  number of repeats within specified windows across 
  the length of each chromosome. This fucntion is designed to work with the merged data output from the GetRepeats function.
}
\usage{
PlotRepeats(merged_csv, window_size = 200000, step_size = 200000)
}
\arguments{
  \item{data}{
    A data frame containing columns \code{Chromosome}, \code{Start}, 
    \code{End}, \code{Match_Start}, and \code{Match_End}, 
    detailing the locations of direct repeats within a genome assembly.
  }
  \item{window_size}{
    An integer specifying the size (in base pairs) of the bin used to count 
    repeats. Defaults to \code{200000}.
  }
  \item{step_size}{
    An integer specifying the step size (in base pairs) between consecutive 
    bins. Defaults to \code{200000}.
  }
}
\value{
  This function does not return an object. Instead, it produces one plot 
  per chromosome, each showing repeat count in windows across 
  the chromosome's length.
}
\author{
  Megan Copeland, Heath Blackmon
}
\references{
  \url{https://github.com/coleoguy/DirectRepeateR}
}
\examples{
# The following code is an example, commented out to prevent execution:
# final_output <- GetRepeats(some_genome_data)
# PlotRepeats(data = merged_csv, window_size = 200000, step_size = 200000)
}

\keyword{visualization}
