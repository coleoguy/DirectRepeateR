\name{PlotRepeats}
\alias{PlotRepeats}
\title{Visualization of Repeat Density Across Chromosomes}
\description{
  The \code{PlotRepeats} function visualizes the density of direct 
  repeats across chromosomes using a sliding window approach. This 
  function generates line plots showing the number of repeats within 
  specified windows across the length of each chromosome.
}
\usage{
PlotRepeats(data, window_size = NULL, step_size = NULL)
}
\arguments{
  \item{data}{
    A data frame generated by the \code{GetRepeats} function, containing 
    details about the locations of direct repeats within a genome assembly.
  }
  \item{window_size}{
    An integer specifying the size of the sliding window used to count 
    the number of repeats. Default is \code{NULL}, which defaults to 
    200,000 base pairs within the function.
  }
  \item{step_size}{
    An integer specifying the step size for the sliding window. This 
    defines the interval at which the window moves along the chromosome. 
    Default is \code{NULL}, which defaults to 200,000 base pairs within 
    the function.
  }
}
\value{
  The function generates and displays line plots for each chromosome, 
  showing the number of repeats within each window along the chromosome's 
  length.
}
\author{
  Megan Copeland, Heath Blackmon
}
\references{
  \url{https://github.com/coleoguy/DirectRepeateR}
}
\examples{
# The following code is an example, but it is commented out to prevent execution.
# PlotRepeats(data = final_output, window_size = 200000, step_size = 200000)
}

\keyword{visualization}
