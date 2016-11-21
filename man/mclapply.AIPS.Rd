\name{mclapply.AIPS}
\alias{mclapply.AIPS}
\title{
Multicore/parallel function that will assign the AIPS models to a gene expression matrix
}
\description{
Given a gene expression matrix D where rows correspond to genes and
columns to samples and a list of Entrez gene ids, this function will
assign all the AIPS models by taking advantage of the multicore/parallel package. 
}
\usage{
mclapply.AIPS (D, EntrezID, n.cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
  A matrix D of gene expression. Rows are genes and columns are samples (tumor from breast cancer patients).
}
  \item{EntrezID}{
  A character vector corresponding to the Entrez Ids of the genes of matrix D.
}
\item{n.cores}{
  Number of CPU (cores) to use to run AIPS. Default is detectCores() - 2
}
}

\details{
  We defined Absolute Inference of Patient Signatures (AIPS) to
  stabilize the assignement of pathway activation and enable assignment
  on single whole-tumor and single cell gene expression. The idea of the
  approach is to use simple rules of the form 
  "if gene A expression < gene B expression assigned sample to subtype X".
  By using those simple rules we are able to assign pathway activation
  using only the expression values of one patient or one single-cell
  profile. We have shown AIPS recapitulates whole-cohort assignments.
  This function will return the pathway activation assignment as well as the posterior probabilities for all the assignments.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cl}{Gene signature activation identified by AIPS. It could be
  either "low", "random/independent", or "high".}
\item{posterior}{A vector corresponding to the posterior probabilities of the
  activation classes in cl.}
\item{gs.info}{A matrix that contains the information about the gene
  signatures (Name, Source, Description)}
\item{raw.aips}{All the individual AIPS assignments}
}
\references{
Paquet ER, Hallett MT.
J Natl Cancer Inst. 2014 Dec 4;107(1):357
}
\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIPS]{aips.models}}, \code{\link[AIPS]{mcgill.gq}}
}
\examples{
## Load the McGill dataset used in the paper
data(mcgill.gq)

## Assignd AIPS on the first 5 gene expression in the McGill dataset
mcgill.AIPS <- mclapply.AIPS (mcgill.gq$D[,1:5],
                              mcgill.gq$EntrezID)

## Print a summary of the assignments for AIPS model #1706
## for all the patients
table(mcgill.AIPS$cl[1706,])

## information of the first model
mcgill.AIPS$gs.info[1706,]

## We can do the samething using only one sample
mcgill.AIPS.first <- mclapply.AIPS (mcgill.gq$D[,1,drop=FALSE], 
                                    mcgill.gq$EntrezID)

## Compare the assignments for patient 1
table(mcgill.AIPS$cl[,1],mcgill.AIPS.first$cl[,1])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
