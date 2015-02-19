\name{read.qPCRfiles}
\alias{read.qPCRfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Read qPCR files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reads in qPCR files for a specific file format.
}
\usage{
read.qPCRfiles(files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
%%     ~~Describe \code{files} here~~
}
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
function (files) 
{
    nfiles <- length(files)
    newdataframe <- NULL
    for (j in 1:nfiles) {
        varnames <- scan(files[[j]], what = "character", nlines = 1, 
            skip = 8, sep = ",", quiet = TRUE)
        tempdata <- scan(files[[j]], what = "character", skip = 9, 
            sep = ",", quiet = TRUE)
        nlines <- floor(length(tempdata)/13)
        datamat <- NULL
        for (i in 1:nlines) {
            temprow <- tempdata[((i - 1) * 13 + 1):(i * 13)]
            datamat <- rbind(datamat, temprow)
        }
        colnames(datamat) <- varnames
        newdatamat <- datamat[datamat[, "Type"] == "Unknown", 
            c("Name", "Type", "Rep. Calc. Conc.")]
        newvals <- as.numeric(newdatamat[newdatamat[, "Rep. Calc. Conc."] != 
            "", "Rep. Calc. Conc."])
        newnames <- as.character(newdatamat[newdatamat[, "Rep. Calc. Conc."] != 
            "", "Name"])
        newdataframe <- rbind(newdataframe, data.frame(Name = newnames, 
            Rep.Calc.Conc = newvals))
    }
    newdataframe[, "Name"] <- as.character(newdataframe[, "Name"])
    return(newdataframe)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
