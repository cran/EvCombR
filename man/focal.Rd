\name{focal}
\alias{focal}
\title{
Focal Elements of a Mass Function
}
\description{
Returns the set of focal elements of a mass function.
}
\usage{
focal(x)
}
\arguments{
  \item{x}{
a mass function
}
}
\value{
focal elements of \code{x}
}
\references{
Dempster, A. P. (1969), A generalization of Bayesian inference, \emph{Journal of the Royal Statistical Society}, \bold{30}, 205-247

Shafer, G., (1976), A Mathematical Theory of Evidence Princeton University Press, 1976
}
\author{
Alexander Karlsson
}

\seealso{
\code{\link{points}}
}
\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass functions 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)

# obtain focal elements
focal(m)
}
