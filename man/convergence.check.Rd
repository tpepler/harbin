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
convergence.check(datavec, graph = c("density", "histogram"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datavec}{
%%     ~~Describe \code{datavec} here~~
}
  \item{graph}{
%%     ~~Describe \code{graph} here~~
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
function (datavec, graph = c("density", "histogram")) 
{
    nbounds <- 7
    boundvals <- quantile(datavec, probs = c(0.2, 0.4, 0.6, 0.8))
    par(mfrow = c(3, 1))
    if (graph[1] == "histogram") {
        hist(datavec, breaks = 20, col = "darkgrey", main = "Gene of interest: Distribution", 
            xlab = "Normalised GOI value", xlim = c(0, max(datavec) + 
                0.2))
    }
    if (graph[1] == "density") {
        d <- density(datavec)
        plot(d, main = "Gene of interest: Distribution", xlab = "Normalised GOI value", 
            xlim = c(0, max(d$x)))
        polygon(d, col = "darkgrey")
    }
    for (i in 1:nbounds) {
        abline(v = boundvals[i], lty = 3)
    }
    nGOIvals <- length(datavec)
    class1UL <- rep(NA, times = nGOIvals - 9)
    class2UL <- rep(NA, times = nGOIvals - 9)
    class3UL <- rep(NA, times = nGOIvals - 9)
    class4UL <- rep(NA, times = nGOIvals - 9)
    class5UL <- rep(NA, times = nGOIvals - 9)
    class1mean <- rep(NA, times = nGOIvals - 9)
    class2mean <- rep(NA, times = nGOIvals - 9)
    class3mean <- rep(NA, times = nGOIvals - 9)
    class4mean <- rep(NA, times = nGOIvals - 9)
    class5mean <- rep(NA, times = nGOIvals - 9)
    for (j in 10:nGOIvals) {
        tempdata <- datavec[1:j]
        tempdata <- tempdata[tempdata >= 0]
        class1UL[j - 9] <- quantile(tempdata, probs = 0.2)
        class2UL[j - 9] <- quantile(tempdata, probs = 0.4)
        class3UL[j - 9] <- quantile(tempdata, probs = 0.6)
        class4UL[j - 9] <- quantile(tempdata, probs = 0.8)
        class5UL[j - 9] <- quantile(tempdata, probs = 1)
        tempboundvals <- quantile(tempdata, probs = c(0, 0.2, 
            0.4, 0.6, 0.8, 1))
        class1data <- tempdata[(tempdata >= tempboundvals[1]) & 
            (tempdata <= tempboundvals[2])]
        class2data <- tempdata[(tempdata > tempboundvals[2]) & 
            (tempdata <= tempboundvals[3])]
        class3data <- tempdata[(tempdata > tempboundvals[3]) & 
            (tempdata <= tempboundvals[4])]
        class4data <- tempdata[(tempdata > tempboundvals[4]) & 
            (tempdata <= tempboundvals[5])]
        class5data <- tempdata[(tempdata > tempboundvals[5]) & 
            (tempdata <= tempboundvals[6])]
        class1mean[j - 9] <- mean(class1data)
        class2mean[j - 9] <- mean(class2data)
        class3mean[j - 9] <- mean(class3data)
        class4mean[j - 9] <- mean(class4data)
        class5mean[j - 9] <- mean(class5data)
    }
    plot(x = (10:nGOIvals), y = class1UL, type = "l", ylim = c(min(rbind(class1UL, 
        class2UL, class3UL, class4UL, class5UL)), max(rbind(class1UL, 
        class2UL, class3UL, class4UL, class5UL))), ylab = "Class upper limits", 
        xlab = "Number of observations", main = "Convergence plot")
    lines(x = (10:nGOIvals), y = class2UL, type = "l")
    lines(x = (10:nGOIvals), y = class3UL, type = "l")
    lines(x = (10:nGOIvals), y = class4UL, type = "l")
    lines(x = (10:nGOIvals), y = class5UL, type = "l")
    plot(x = (10:nGOIvals), y = class1mean, type = "l", ylim = c(min(rbind(class1mean, 
        class2mean, class3mean, class4mean, class5mean)), max(rbind(class1mean, 
        class2mean, class3mean, class4mean, class5mean))), ylab = "Class means", 
        xlab = "Number of observations", main = "Convergence plot")
    lines(x = (10:nGOIvals), y = class2mean, type = "l")
    lines(x = (10:nGOIvals), y = class3mean, type = "l")
    lines(x = (10:nGOIvals), y = class4mean, type = "l")
    lines(x = (10:nGOIvals), y = class5mean, type = "l")
    par(mfrow = c(1, 1))
    return(boundvals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line