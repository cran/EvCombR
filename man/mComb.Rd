\name{mComb}
\alias{mComb}
\title{
Modified Dempster's Combination Operator
}
\description{
Combine evidence in the form of mass functions using modified Dempster's combination operator.
}
\usage{
mComb(x,y,z)
}

\arguments{
\item{x}{
    single mass function or a list of mass functions    
}
\item{y}{
    single mass function if \code{x} is a single mass function, a prior distribution  or missing if \code{x} is a list 
}
\item{z}{
    prior distribution if \code{x} and \code{y} are mass functions, otherwise missing
}
}

\details{
The prior distribution is provided in the form of a list where the names are equivalent to the state space. See the examples.
}

\value{
mass function
}

\references{
Fixsen, D., Mahler, R. P. S. (1997), The modified Dempster-Shafer approach to classification, \emph{IEEE Transactions on Systems, Man and Cybernetics, Part A: Systems and Humans}, \bold{27}, 96-104
}

\author{
Alexander Karlsson
}


\seealso{
\code{\link{dComb}, \link{yComb}, \link{cComb}}
}   

\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass functions 
m1 <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)
m2 <- mass(list("a"=0.2, "a/b/c"=0.8), stateSpace)

# modified Dempster's combination using the uniform prior
mComb(m1, m2)
# or 
mComb(list(m1, m2))     

# modified Dempster's combination using a specific prior
mComb(m1, m2, list("a"=0.1, "b"=0.1, "c"=0.8))
# or 
mComb(list(m1, m2), list("a"=0.1, "b"=0.1, "c"=0.8)) 
}

