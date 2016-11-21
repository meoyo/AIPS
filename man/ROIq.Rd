\name{ROIq}
\alias{ROIq}
\title{
Function that will assign the AIPS models to a gene expression matrix
}
\description{
Given a gene expression matrix D where rows correspond to genes and
columns to samples, and a gene set, this function will compute the ROIq
assignements (low, independent, or high) associated to samples in matrix
D.
}
\usage{
ROIq (D,genename,gene.set,n=1000,q=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
  A matrix D of gene expression. Rows are genes and columns are samples (tumor from breast cancer patients).
}
  \item{genename}{
  A character vector corresponding to the genenames of the genes of
  matrix D. Expect length(genename) == nrow(D)
}
\item{gene.set}{ A list element (list(up=,down=)) corresponding to the
  genes in the gene signature. The up-regulated genes should be set in
  "up" and the down-regulated genes should be set in "down".
}
\item{n}{ Number of random sampling necessary to define the region of independence.
}
\item{q}{ The fraction of random sampling contained within the region of
independence.} 
}

\details{
  This function will compute the Region Of Independence at quantile q
  (ROIq) descibe in Paquet et al. Given a gene expression dataset, a
  list of gene names, and a gene signature the function will identify the
  region of independence and will return assignments associated to the
  patients given the gene signatures. The assignment could be low,
  independent, or high in function of how strongly the patients are
  associated to the gene signature.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cl}{Gene signature activation identified by AIPS. It could be
  either "low", "random/independent", or "high".}
}
\references{
Paquet ER, Hallett MT.
J Natl Cancer Inst. 2014 Dec 4;107(1):357
}
\author{
  Eric R. Paquet (eric.r.paquet@gmail.com),
  Ali Tofigh (alix.tofigh@gmail.com),
  Robert Lesurf (robert.lesurf@gmail.com)
}

\seealso{
\code{\link[AIPS]{mcgill.gq}}
}
\examples{
## Load the McGill dataset used in the paper
data(mcgill.gq)
data(aips.models)

## Get one esr1 (estrogen) gene signature
esr1.gs = aips.models[[1706]]$gs.bresat

## Assigned the ROI95 on McGill dataset using the esr1 gene signature
mcgill.ROI95.esr1.gs <- ROIq (mcgill.gq$D,
                              mcgill.gq$EntrezID,
                              list(up=esr1.gs$ENTREZ$up,down=esr1.gs$ENTREZ$down))

## Print a summary of all the assignments
table(mcgill.ROI95.esr1.gs$cl)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
