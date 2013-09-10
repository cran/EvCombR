\name{dComb}
\alias{dComb}

\title{
Dempster's Combination Operator
}
\description{
Combine evidence in the form of mass functions using Dempster's combination operator. 
}
\usage{
dComb(x,y)
}

\arguments{
\item{x}{
    single mass function or a list of mass functions
}
\item{y}{
    single mass function if \code{x} is a single mass function, otherwise missing    

}
}

\value{
mass function
}

\references{
Dempster, A. P. (1969), A generalization of Bayesian inference, 
\emph{Journal of the Royal Statistical Society}, \bold{30}, 205-247

Shafer, G. (1976), A Mathematical Theory of Evidence Princeton University Press
}

\author{
Alexander Karlsson
}

\seealso{
\code{\link{yComb}, \link{mComb}, \link{cComb}}
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass functions 
m1 <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)
m2 <- mass(list("a"=0.2, "a/b/c"=0.8), stateSpace)

# Dempster's combination
dComb(m1, m2)
# or 
dComb(list(m1, m2))
}

