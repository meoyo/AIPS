\name{aips.models}
\alias{aips.models}
\title{
All the AIPS models.
}
\description{
This is the data object that contains all the AIPS models. It contains
the naive bayes classifier using the selected pairs of genes for all the
trained models presented in Paquet et al.
}
\usage{
aips.models
}
%- maybe also 'usage' for other objects documented here.

\details{
  This is the data object that contains all the AIPS models.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{all.pairs}{The rules used in the model definition form EntrezID gene A < EntrezID gene B}
\item{k}{The selected number of optimal rules.}
\item{one.vs.all.tsp}{The Naive bayes classifier used in combination
  with the all.pairs rules}
\item{selected.pairs.list}{The list of rules sorted from the best discriminating rule to the least discriminating rules subdivided by assignments.}
}
\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIPS]{apply.AIPS}}, \code{\link[AIPS]{mcgill.gq}}, \code{\link[e1071]{naiveBayes}}
}

\examples{
## Load a sample of the McGill dataset used in the paper
data(aips.models)

## list the number of models
length(aips.models)

## list the model
aips.models[[1]]$model$model

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
