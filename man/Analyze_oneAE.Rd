\name{Analyze_oneAE}
\alias{Analyze_oneAE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model selection in logistic regression using Hastings-Metropolis algorithm
}
\description{
In this function, we propose to use a logistic regression whose sparsity is viewed as a model selection
challenge. Since the model space is huge, a Metropolis-Hastings algorithm carries out the model selection
by maximizing the BIC criterion
}
\usage{
Analyze_oneAE(ae, drug, maxit, alpha, nbinit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ae}{Binary vector indicate if individual suffers from adverse event (1) or no (0). 
%%     ~~Describe \code{ae} here~~
}
  \item{drug}{Matrix of drugs consumptions. Each row corresponds to one individual drug consumptions. 
%%     ~~Describe \code{drug} here~~
}
  \item{maxit}{
  Numeric indicating the number of iterations.
%%     ~~Describe \code{maxit} here~~
}
  \item{alpha}{
  Numeric indicating the size of the neighborhood where the proposal will be
  uniformly sampled at each iteration. See~\url{http://arxiv.org/abs/1505.03366} 
  for more details.
%%     ~~Describe \code{alpha} here~~
}
  \item{nbinit}{
  The number of random initializations of the algorithm.
%%     ~~Describe \code{nbinit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of detected signals.
}
\references{
Matthieu Marbac, Pascale Tubert-Bitter, Mohammed Sedki

Bayesian model selection in logistic regression for the detection of adverse drug reactions~\url{http://arxiv.org/abs/1505.03366}
}
\author{
Matthieu Marbac and Mohammed Sedki
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Binary data} 
\keyword{Logistic regression}
\keyword{Metropolis-Hastings}
\keyword{Model selection} 
\keyword{Pharmacovigilance} 
\keyword{Spontaneous reporting}

