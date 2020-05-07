\name{ATAC_JV_adipose}
\alias{ATAC_JV_adipose}
\docType{data}
\title{
    ATAC JV adipose data
}
\description{
  A data set for testing the iregnet package. Specifically features of neuroblastoma
  patients and their corresponding survival times in interval format.
}
\usage{data("ATAC_JV_adipose")}
\format{
  The format is:
  
  Inputs are a input matrix of 465 rows/observations and 29 columns/features.
  The outputs are in order however row names are applied to match inputs and
  outputs, just in case. The form of outputs are 465 rows/observations and 2
  columns/surv time, which is an interval of survival time for each observation
}
\details{
  data( ATAC_JV_adipose , package="iregnet")
}
\source{
https://github.com/andruuhurst/neuroblastoma-data/tree/master/data/ATAC_JV_adipose
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ATAC_JV_adipose)
## maybe str(ATAC_JV_adipose) ; plot(ATAC_JV_adipose) ...
}
\keyword{datasets}
