\name{PlotRepeats}
\alias{PlotRepeats}
\title{Visualization of Repeat Counts Across Chromosomes}
\description{
  The \code{PlotRepeats} function visualizes the distribution of direct 
  repeats across chromosomes by generating line plots 
  showing the number of repeats within specified windows across 
  the length of each chromosome. This function is designed to work with 
  the merged data output from the \code{GetRepeats} function.
}
\usage{
PlotRepeats(data, window_size = 200000, step_size = 200000)
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
  Invisibly returns a named list of \pkg{ggplot2} plot objects, one per 
  chromosome, so they can be saved or modified. The plots are also 
  printed, one per chromosome, showing repeat counts in windows across 
  each chromosome's length.
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
# plots <- PlotRepeats(data = final_output, window_size = 200000, step_size = 200000)
}

\keyword{visualization}
