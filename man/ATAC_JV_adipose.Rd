\name{ATAC_JV_adipose}
\alias{ATAC_JV_adipose}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("ATAC_JV_adipose")}
\format{
  The format is:
List of 2
 $ inputs : num [1:465, 1:29] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:465] "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr10:18024675-38818835" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr11:69774695-87688378" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr12:145739-7189876" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr14:19050000-107289540" ...
  .. ..$ : chr [1:29] "quartile.0%" "quartile.25%" "quartile.50%" "quartile.75%" ...
 $ outputs: num [1:465, 1:2] 8.59 7.84 10.08 8.19 6.87 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:465] "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr10:18024675-38818835" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr11:69774695-87688378" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr12:145739-7189876" "ATAC_JV_adipose/samples/AC1/MSC77/problems/chr14:19050000-107289540" ...
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
data(ATAC_JV_adipose)
## maybe str(ATAC_JV_adipose) ; plot(ATAC_JV_adipose) ...
}
\keyword{datasets}
