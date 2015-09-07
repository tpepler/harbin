\name{harbin.test}
\alias{harbin.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Harbin test (EXPERIMENTAL)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs the Harbin test for two groups. (PLEASE NOTE that this method is experimental and still undergoing testing.)
}
\usage{
harbin.test(x, y, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{First dataset, as a vector.}
  \item{y }{Second dataset, as a vector.}
  \item{reps }{Number of bootstrap replications to use.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The Harbin test is a non-parametric test for two-sample location-scale-shape problem, testing the hypothesis that two datasets originated from population distributions which can be described by the same probability distribution function. The alternative hypothesis is that the location, variability and/or shapes of the two population distributions differ.

The test statistic is calculated from quantiles of the pooled datasets, and compared to a bootstrap distribution of the statistic under the null hypothesis.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the following components:
\item{statistic }{Test statistic (the proportion of samples in the first dataset for which the labels have changed when the second dataset was added.}
\item{crit.val }{95\% critical value for the test.}
\item{p.value }{P-value for the test.}
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
# Example 1: Difference in location
data1 <- rnorm(n = 30, mean = 0, sd = 1)
data2 <- rnorm(n = 20, mean = 0.5, sd = 1)
harbin.test(x = data1, y = data2)

# Example 2: Difference in location and variability
data1 <- rnorm(n = 30, mean = 0, sd = 1)
data2 <- rnorm(n = 20, mean = 0.2, sd = 1.5)
harbin.test(x = data1, y = data2)

# Example 3: Difference in location and shape
data1 <- rnorm(n = 30, mean = 0, sd = 1)
data2 <- runif(n = 20, min = -2.5, max = 3.5)
harbin.test(x = data1, y = data2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
