\name{ConvertToGFF}
\alias{ConvertToGFF}
\title{Convert Repeat Detection Output to GFF Format}
\description{
  The \code{ConvertToGFF} function converts the output from a repeat detection process into a GFF (General Feature Format) file. It processes each row of the input data frame, representing detected repeat sequences, and generates GFF entries for the full repeat region, the first copy of the repeat, and the second copy of the repeat.
}
\usage{
  ConvertToGFF(data)
}
\arguments{
  \item{data}{
    A data frame containing the output from the repeat detection process. Each row should represent a detected repeat, with columns for the chromosome, start position, end position, match start position, and match end position.
  }
}
\value{
  The function writes a GFF file to the fixed location \code{"results/gff_output.txt"}. It does not return a value.
}
\author{
  Megan Copeland, Heath Blackmon
}
\examples{
# The following code is an example, but it is commented out to prevent execution.
# ConvertToGFF(data = final_output)
}
\keyword{data manipulation}