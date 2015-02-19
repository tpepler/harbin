\name{GOI.normalise}
\alias{GOI.normalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Gene of interest normalisation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Normalises gene of interest data using selected reference genes.
}
\usage{
GOI.normalise(GOIfiles, refgenefiles, refdatabase = NULL, write.output = FALSE, pathname = "", filename = "GOI_output.csv", graph = c("density", "histogram"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOIfiles}{
%%     ~~Describe \code{GOIfiles} here~~
List containing the names (addresses) of the gene of interest files.
}
  \item{refgenefiles}{
%%     ~~Describe \code{refgenefiles} here~~
List of lists containing the names (addresses) of the reference gene files.
}
  \item{refdatabase}{
%%     ~~Describe \code{refdatabase} here~~
Address of the reference database file (if applicable).
}
  \item{write.output}{
%%     ~~Describe \code{write.output} here~~
Logical (TRUE/FALSE), indicating whether output should be written to a reference database file.
}
  \item{pathname}{
%%     ~~Describe \code{pathname} here~~
System path of the reference data base file (if refdatabase is not NULL).
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
Name of the reference data base file to be written out (applicable if write.output = TRUE).
}
  \item{graph}{
%%     ~~Describe \code{graph} here~~
Type of graph to construct: Density plot (default) or a histogram.
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
P.T. Pepler (tpepler@sun.ac.za)
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
function (GOIfiles, refgenefiles, refdatabase = NULL, write.output = FALSE, 
    pathname = "", filename = "GOI_output.csv", graph = c("density", 
        "histogram")) 
{
    nrefgenes <- length(refgenefiles)
    GOIdata <- read.qPCRfiles(GOIfiles)
    n <- nrow(GOIdata)
    refgenedata <- vector("list", nrefgenes)
    refgenevals <- vector("list", nrefgenes)
    for (j in 1:nrefgenes) {
        refgenedata[[j]] <- read.qPCRfiles(refgenefiles[[j]])
        for (i in 1:n) {
            if (GOIdata[i, "Name"] != refgenedata[[j]][i, "Name"]) {
                stop("Warning: Reference gene names do not match gene of interest names!\n")
            }
        }
        refgenevals[[j]] <- refgenedata[[j]][, "Rep.Calc.Conc"]
    }
    refindex <- (Reduce("*", refgenevals))^(1/nrefgenes)
    GOI.normdata <- GOIdata[, "Rep.Calc.Conc"]/refindex
    refbaseadd.answer <- "n"
    if (!is.null(refdatabase)) {
        refdatabase.data <- read.csv(paste(pathname, refdatabase, 
            sep = ""), head = TRUE)
        for (i in 1:n) {
            if (GOIdata[i, "Name"] \%in\% refdatabase.data[, "Name"]) {
                cat(paste("Warning: GOI name ", GOIdata[i, "Name"], 
                  " already found in reference database!\n", 
                  sep = ""))
            }
        }
        boundvals <- convergence.check(c(refdatabase.data[, "GOI.normalised"], 
            GOI.normdata), graph = graph)
        refdatabase.index <- rep(0, times = nrow(refdatabase.data))
        refdatabase.index[((refdatabase.data[, "GOI.normalised"] >= 
            0) & (refdatabase.data[, "GOI.normalised"] <= boundvals[1]))] <- 1
        refdatabase.index[((refdatabase.data[, "GOI.normalised"] > 
            boundvals[1]) & (refdatabase.data[, "GOI.normalised"] <= 
            boundvals[2]))] <- 2
        refdatabase.index[((refdatabase.data[, "GOI.normalised"] > 
            boundvals[2]) & (refdatabase.data[, "GOI.normalised"] <= 
            boundvals[3]))] <- 3
        refdatabase.index[((refdatabase.data[, "GOI.normalised"] > 
            boundvals[3]) & (refdatabase.data[, "GOI.normalised"] <= 
            boundvals[4]))] <- 4
        refdatabase.index[refdatabase.data[, "GOI.normalised"] > 
            boundvals[4]] <- 5
        refdatabase.data[, "Interval"] <- refdatabase.index
        harbin.out <- harbin.test(x = refdatabase.data[, "GOI.normalised"], 
            y = GOI.normdata)
        cat(paste("\nProportion of labels changing in reference data base: ", 
            round(harbin.out$statistic * 100, 1), "\%\n", sep = ""))
        cat("H0: New data originated from same distribution as reference data\n")
        cat("H1: New data and reference data come from different distributions\n")
        cat(paste("Harbin test p-value = ", harbin.out$p.value, 
            "\n", sep = ""))
        if (harbin.out$p.value <= 0.05) {
            cat("\nWARNING: Reference data and new data may not be compatible!\n")
        }
        cat("\nDo you want to add the new data to the reference database (y/n): ")
        refbaseadd.answer <- readline()
        cat("\n")
    }
    else {
        boundvals <- convergence.check(GOI.normdata, graph = graph)
    }
    interval.index <- rep(0, times = n)
    interval.index[((GOI.normdata >= 0) & (GOI.normdata <= boundvals[1]))] <- 1
    interval.index[((GOI.normdata > boundvals[1]) & (GOI.normdata <= 
        boundvals[2]))] <- 2
    interval.index[((GOI.normdata > boundvals[2]) & (GOI.normdata <= 
        boundvals[3]))] <- 3
    interval.index[((GOI.normdata > boundvals[3]) & (GOI.normdata <= 
        boundvals[4]))] <- 4
    interval.index[GOI.normdata > boundvals[4]] <- 5
    resultsmat <- data.frame(Name = GOIdata[, "Name"], GOI.data = GOIdata[, 
        "Rep.Calc.Conc"], Reference.index = refindex, GOI.normalised = GOIdata[, 
        "Rep.Calc.Conc"]/refindex, Interval = interval.index)
    if (write.output) {
        write.csv(resultsmat, paste(pathname, filename, sep = ""), 
            row.names = FALSE, quote = FALSE)
    }
    if ((!is.null(refdatabase)) & (refbaseadd.answer \%in\% c("y", 
        "Y", "yes", "Yes", "YES"))) {
        write.csv(rbind(refdatabase.data, resultsmat), paste(pathname, 
            refdatabase, sep = ""), row.names = FALSE, quote = FALSE)
    }
    if (is.null(refdatabase)) {
        cat("Gene of interest results:\n\n")
    }
    else {
        cat("Gene of interest (new data) results:\n\n")
    }
    return(resultsmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
