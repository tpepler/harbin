\name{open.files}
\alias{open.files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Open files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Graphical user interface for selecting files to open.
}
\usage{
open.files()
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
    fileName <- tclvalue(tkgetOpenFile(multiple = T))
    if (!nchar(fileName)) {
        tkmessageBox(message = "No file was selected!")
        return(NULL)
    }
    endpoints <- str_locate_all(fileName, "}")[[1]]
    nfiles <- nrow(endpoints)
    filelist <- vector("list", nfiles)
    for (i in 1:nfiles) {
        if (i == 1) {
            filelist[i] <- substr(fileName, 2, endpoints[i, "end"] - 
                1)
        }
        if (i > 1) {
            filelist[i] <- substr(fileName, endpoints[(i - 1), 
                "end"] + 3, endpoints[i, "end"] - 1)
        }
    }
    return(filelist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
