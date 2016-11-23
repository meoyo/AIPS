\name{aips.normal.info}
\alias{aips.normal.info}
\title{
Information about AIPS assignments of normal-like samples in the
validation set (n = 1261).
}
\description{
This is the data object that contains all the information about the
normal-like assignments in the validation set (n = 1261). Using this
data.frame it is possible to compare your samples' assignments to those
of normal-like samples.
}
\usage{
aips.normal.info
}
%- maybe also 'usage' for other objects documented here.

\details{
  This is the data object that contains information about the AIPS
  assignments for Normal-like samples in the validation set (n=1261).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Name}{Name of the AIPS model}
\item{Source}{Source of the gene signature used to build the AIPS modes}
\item{Description}{More information on the AIPS model}
\item{Normal.Assignment}{The main assignment for the Normal-like
  samples}
\item{Normal.low.percentage}{Percentage of Normal-like patients assigned
  to the low class}
\item{Normal.random.percentage}{Percentage of Normal-like patients assigned
  to the random/independent class}
\item{Normal.high.percentage}{Percentage of Normal-like patients assigned
  to the high class}
}
\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIPS]{apply.AIPS}}, \code{\link[AIPS]{mcgill.gq}}
}

\examples{
## Load a sample of the McGill dataset used in the paper
data(aips.normal.info)

## Show first AIPS models information
head(aips.normal.info)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
