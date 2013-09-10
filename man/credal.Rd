\name{credal}
\alias{credal}

\title{
Constructor Function for Credal Sets (based on intervals)
}
\description{
Construct a credal set based on probability intervals or a single probability function. The algorithm used for finding the 
extreme points corresponding to lower and upper bounds is described in De Campos et al. (1994).
}
\usage{
credal(x, y, z)
}

\arguments{
\item{x}{
lower bounds of probabilities or a single probability function (in the form of a numeric vector)
}
\item{y}{
character vector representing the state space if \code{x} is a single probability function or upper bounds for probabilities if \code{x} are lower bounds 
}
\item{z}{
character vector representing the state space if \code{x} and \code{y} are lower respective upper bounds of probabilities
}
}

\value{
A credal set represented by a set of extreme points.
}
\references{
Levi, I. (1983), The enterprise of knowledge, The MIT press

Arnborg, S. (2006), Robust Bayesianism: Relation to Evidence Theory, \emph{Journal of Advances in Information Fusion}, \bold{1}, 63-74

Karlsson, A., Johansson, R., Andler, S. F. (2011), Characterization and Empirical Evaluation of Bayesian and Credal Combination Operators, \emph{Journal of Advances in Information Fusion}, \bold{6}, 150-166

De Campos L. M., Huete, J. F., Moral S., Probability Intervals: a Tool for Uncertain Reasoning,\emph{International Journal of Uncertainty, Fuzziness, and Knowledge-Based Systems}, \bold{2}, 167-196
}

\author{
Alexander Karlsson
}


\seealso{
\code{\link{cComb}}
}
\examples{
# state space
stateSpace <- c("a", "b", "c")

# lower and upper bounds
c1 <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace)

# single probability function
c2 <- credal(c(0.1, 0.2, 0.7), stateSpace)

}

