\name{H3K36me3_TDH_immune}
\alias{H3K36me3_TDH_immune}
\docType{data}
\title{
   H3K36me3 TDH immune data
}
\description{
  A data set for testing the iregnet package. Specifically features of neuroblastoma
  patients and their corresponding survival times in interval format.
}
\usage{data("H3K36me3_TDH_immune")}
\format{
  The format is:

  Inputs are a input matrix of 63 rows/observations and 28 columns/features.
  The outputs are in order however row names are applied to match inputs and
  outputs, just in case. The form of outputs are 63 rows/observations and 2
  columns/surv time, which is an interval of survival time for each observation
}
\details{
    data( H3K36me3_TDH_immune , package="iregnet")
}
\source{
    https://github.com/andruuhurst/neuroblastoma-data/tree/master/data/H3K36me3_TDH_immune
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(H3K36me3_TDH_immune)
## maybe str(H3K36me3_TDH_immune) ; plot(H3K36me3_TDH_immune) ...
}
\keyword{datasets}
