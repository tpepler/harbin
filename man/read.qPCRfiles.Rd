\name{read.qPCRfiles}
\alias{read.qPCRfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read RT-qPCR files (Qiagen Rotor-Gene Q thermal cycler)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reads in RT-qPCR files produced by the Qiagen Rotor-Gene Q thermal cycler from the local file system.
}
\usage{
read.qPCRfiles(files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files }{List containing the local file system address(es) of RT-qPCR file(s) to be read.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a utility function used by \code{\link{harbin.quant}} as part of the graphical user interface to perform Harbin quantification.
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