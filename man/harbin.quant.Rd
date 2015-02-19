\name{harbin.quant}
\alias{harbin.quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Harbin quantification
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs Harbin quantification for RT-qPCR data.
}
\usage{
harbin.quant()
}
%- maybe also 'usage' for other objects documented here.
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    require(tcltk)
    require(stringr)
    tkmessageBox(message = "Select the gene of interest file(s)")
    GOIfiles <- open.files()
    refgenefiles <- vector("list", 3)
    tkmessageBox(message = "Select the Reference gene 1 file(s)")
    refgenefiles[[1]] <- open.files()
    tkmessageBox(message = "Select the Reference gene 2 file(s)")
    refgenefiles[[2]] <- open.files()
    tkmessageBox(message = "Select the Reference gene 3 file(s)")
    refgenefiles[[3]] <- open.files()
    GOI.normalise(GOIfiles, refgenefiles)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
