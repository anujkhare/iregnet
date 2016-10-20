\name{penalty.learning}
\alias{penalty.learning}
\docType{data}
\title{
  Penalty learning data
}
\description{
  A data set for testing the iregnet package. The problem is learning
  a penalty function for predicting the optimal number of peaks in
  ChIP-seq data.
}
\usage{data("penalty.learning")}
\format{
X.mat is a numeric matrix with 443 rows/observations and 22
columns/features. Each observation is a genomic segmentation problem,
for which the features (quartiles, mean, etc) are easily computed. The
y.mat is a numeric matrix with 2 columns. The first column is the
minimal log penalty value that predicts a minimal error model (and
-Inf if there is no lower limit). The second column is the maximum log
penalty value that predicts a minimal error model (and Inf if there is
no upper limit). The goal is to learn a penalty function f such that
lower.limit < f(X.mat) < upper.limit.
}
\source{
  data(H3K4me3.PGP.immune.4608, package="PeakSegJoint")
}
\references{
  http://jmlr.org/proceedings/papers/v28/hocking13.html
}
\keyword{datasets}
