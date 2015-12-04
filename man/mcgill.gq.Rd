\name{mcgill.gq}
\alias{mcgill.gq}
\title{
Sample of a breast cancer gene expression dataset generated at McGill University
}
\description{
This is a breast cancer gene expression dataset generated at McGill University using the Affymetrix Gene ST platform. This is the validation dataset used in the A-PAM50 paper and in use in the examples of AIPS.
}
\usage{
mcgill.gq
}
%- maybe also 'usage' for other objects documented here.

\details{
A sample of a breast cancer gene expression dataset currently generated at McGill University on the Affymetrix gene ST platform.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{D}{A matrix of gene expression values, rows are genes and columns are samples}
\item{EntrezID}{A vector of Entrez ids corresponding to genes (rows) of D}
}

\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIPS]{aips.models}}, \code{\link[AIPS]{apply.AIPS}}
}
\examples{
## Load a sample of the McGill dataset used in the paper
data(mcgill.gq)

## Dimensions of the gene expression matrix 
dim(mcgill.gq$D)

## Number of EntrezID
length(mcgill.gq$EntrezID)
}
\keyword{ datasets }
