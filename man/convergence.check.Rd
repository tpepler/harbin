\name{convergence.check}
\alias{convergence.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Convergence check
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to check if the Harbin quantification boundaries converge.
}
\usage{
%convergence.check(datavec, graph = c("density", "histogram"))
convergence.check(datavec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datavec }{Vector containing normalised gene-of-interest (GOI) values.}
  %\item{graph }{Type of graph ("density", "histogram") to be used for represenation of distribution of normalised GOI values (default = "density").}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Gives a graphical estimator of the distribution of the GOI values, together with plots detailing the change in the cut-off values for Harbin quantification as each additional sample is added to the dataset. Sufficient convergence of these cut-off values should happen at the sample size which is sufficient to describe the distribution of the GOI values accurately.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Theo Pepler
}
\note{
%%  ~~further notes~~
This is a utility function used by \code{\link{harbin.quant}}, but can also be called directly if required.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{harbin.quant}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
