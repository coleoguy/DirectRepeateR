\name{PlotRepeats}
\alias{PlotRepeats}
\title{Visualization of Rolling-Average Repeat Counts Across Chromosomes}
\description{
  The \code{PlotRepeats} function visualizes the distribution of direct 
  repeats across chromosomes, using a sliding window and applying a rolling 
  average to smooth out fluctuations. This function generates line plots 
  showing the smoothed number of repeats within specified windows across 
  the length of each chromosome. This fucntion is designed to work with the merged data output from the GetRepeats function.
}
\usage{
PlotRepeats(merged_csv, window_size = 200000, step_size = 200000, rolling_window = 5)
}
\arguments{
  \item{data}{
    A data frame containing columns \code{Chromosome}, \code{Start_Position}, 
    \code{End_Position}, \code{Match_Position}, and \code{Match_End_Position}, 
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
  \item{rolling_window}{
    An integer specifying the number of bins over which to compute a rolling 
    average. Defaults to \code{5}, meaning each bin is averaged with the 2 
    bins on either side, if available.
  }
}
\value{
  This function does not return an object. Instead, it produces one plot 
  per chromosome, each showing the rolling-average repeat count across 
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
# PlotRepeats(data = merged_csv, window_size = 200000, step_size = 200000, rolling_window = 5)
}

\keyword{visualization}
