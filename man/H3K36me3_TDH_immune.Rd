\name{H3K36me3_TDH_immune}
\alias{H3K36me3_TDH_immune}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("H3K36me3_TDH_immune")}
\format{
  The format is:
List of 2
 $ inputs : num [1:63, 1:28] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:63] "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr10:18024675-38818835" "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr10:51448845-125869472" "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr12:7239876-34856694" "H3K36me3_TDH_immune/samples/monocyte/McGill0001/problems/chr10:18024675-38818835" ...
  .. ..$ : chr [1:28] "quartile.0%" "quartile.25%" "quartile.50%" "quartile.75%" ...
 $ outputs: num [1:63, 1:2] 10.9 10.4 9.3 11.9 10.5 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:63] "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr10:18024675-38818835" "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr10:51448845-125869472" "H3K36me3_TDH_immune/samples/bcell/McGill0322/problems/chr12:7239876-34856694" "H3K36me3_TDH_immune/samples/monocyte/McGill0001/problems/chr10:18024675-38818835" ...
  .. ..$ : chr [1:2] "min.log.lambda" "max.log.lambda"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(H3K36me3_TDH_immune)
## maybe str(H3K36me3_TDH_immune) ; plot(H3K36me3_TDH_immune) ...
}
\keyword{datasets}
