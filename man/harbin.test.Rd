\name{harbin.test}
\alias{harbin.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Harbin test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs the Harbin test for two groups.
}
\usage{
harbin.test(x, y, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
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
function (x, y, reps = 1000) 
{
    x.n <- length(x)
    y.n <- length(y)
    bound <- quantile(x, probs = c(0.2, 0.4, 0.6, 0.8))
    indvec <- rep(1, times = length(x))
    indvec[((x > bound[1]) & (x <= bound[2]))] <- 2
    indvec[((x > bound[2]) & (x <= bound[3]))] <- 3
    indvec[((x > bound[3]) & (x <= bound[4]))] <- 4
    indvec[x > bound[4]] <- 5
    changeprop <- rep(NA, times = reps)
    for (i in 1:reps) {
        bootsamp <- x[sample(1:x.n, size = y.n, replace = TRUE)]
        bootbound <- quantile(c(x, bootsamp), probs = c(0.2, 
            0.4, 0.6, 0.8))
        bootindvec <- rep(1, times = x.n)
        bootindvec[((x > bootbound[1]) & (x <= bootbound[2]))] <- 2
        bootindvec[((x > bootbound[2]) & (x <= bootbound[3]))] <- 3
        bootindvec[((x > bootbound[3]) & (x <= bootbound[4]))] <- 4
        bootindvec[x > bootbound[4]] <- 5
        changevec <- bootindvec - indvec
        changeprop[i] <- (x.n - length(changevec[changevec == 
            0]))/x.n
    }
    newbound <- quantile(c(x, y), probs = c(0.2, 0.4, 0.6, 0.8))
    newindvec <- rep(1, times = x.n)
    newindvec[((x > newbound[1]) & (x <= newbound[2]))] <- 2
    newindvec[((x > newbound[2]) & (x <= newbound[3]))] <- 3
    newindvec[((x > newbound[3]) & (x <= newbound[4]))] <- 4
    newindvec[x > newbound[4]] <- 5
    newchangevec <- newindvec - indvec
    newchangeprop <- (x.n - length(newchangevec[newchangevec == 
        0]))/x.n
    p.value <- length(changeprop[changeprop >= newchangeprop])/reps
    return(list(statistic = newchangeprop, crit.val = quantile(changeprop, 
        prob = 0.95), p.value = p.value))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
