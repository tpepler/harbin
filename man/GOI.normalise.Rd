\name{GOI.normalise}
\alias{GOI.normalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Gene of interest normalisation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Utility function to normalise gene of interest data using selected reference genes.
}
\usage{
%GOI.normalise(GOIfiles, refgenefiles, refdatabase = NULL, write.output = FALSE, pathname = "", filename = "GOI_output.csv", graph = c("density", "histogram"))
GOI.normalise(GOIfiles, refgenefiles, refdatabase = NULL, write.output = FALSE, pathname = "", filename = "GOI_output.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOIfiles }{List containing the names (addresses) of the gene of interest files.}
  \item{refgenefiles }{List of lists containing the names (addresses) of the reference gene files.}
  \item{refdatabase }{Address of the reference database file (if applicable).}
  \item{write.output }{Logical (TRUE/FALSE), indicating whether output should be written to a reference database file.}
  \item{pathname }{System path of the reference data base file (if refdatabase is not NULL).}
  \item{filename }{Name of the reference data base file to be written out (applicable if write.output = TRUE).}
  %\item{graph }{Type of graph to construct: Density plot (default) or a histogram.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
