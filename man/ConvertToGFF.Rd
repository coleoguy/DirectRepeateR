\name{ConvertToGFF}
\alias{ConvertToGFF}
\title{Convert Repeat Detection Output to GFF3 Format}
\description{
  The \code{ConvertToGFF} function converts the output from the 
  \code{GetRepeats} function into GFF3 (General Feature Format version 3). 
  For each detected repeat it generates three entries: a 
  \code{repeat_region} spanning the full repeat (start of the first copy 
  to the end of the second copy) and a \code{direct_repeat} feature for 
  each of the two copies, linked to the region via \code{Parent} 
  attributes.
}
\usage{
  ConvertToGFF(data, path = NULL)
}
\arguments{
  \item{data}{
    A data frame containing the output from the repeat detection process. 
    Each row should represent a detected repeat, with columns 
    \code{Chromosome}, \code{Start}, \code{End}, \code{Match_Start}, and 
    \code{Match_End}.
  }
  \item{path}{
    (Optional) A character string specifying the path for a GFF3 output 
    file. If \code{NULL} (the default), no file is written and the GFF 
    table is only returned. The written file starts with the 
    \code{##gff-version 3} pragma and contains no column header row, so 
    it can be read by IGV, rtracklayer, bedtools, and similar tools.
  }
}
\value{
  A nine-column data frame in GFF3 column layout 
  (\code{seqid}, \code{source}, \code{type}, \code{start}, \code{end}, 
  \code{score}, \code{strand}, \code{phase}, \code{attributes}).
}
\author{
  Megan Copeland, Heath Blackmon
}
\examples{
# The following code is an example, but it is commented out to prevent execution.
# gff <- ConvertToGFF(data = final_output)
# gff <- ConvertToGFF(data = final_output, path = "repeats.gff3")
}
\keyword{data manipulation}
