\name{mass}
\alias{mass}
\title{
Constructor Function for Mass Functions
}
\description{
Construct a mass function based on a named list of focal elements or a \code{\link{massQ-class}} object. For more information, see the details section. 
}
\usage{
mass(x, y)
}

\arguments{
\item{x}{
a named list of focal elements or a \code{\link{massQ-class}} object            
}
\item{y}{
a character vector representing the state space or missing if \code{x} is an \code{massQ} object. 
}
}
\details{
Focal elements are represented by the notation "<\emph{s1}>/.../<\emph{sn}>" where <\emph{s1}>...<\emph{sn}> are any states within the state space (see the examples below). Note that the word "ES" and the symbol "/" are reserved.
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
\code{\link{dComb}, \link{mComb}, \link{yComb}}
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# construct mass functions 
m1 <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)
m2 <- mass(list("a"=0.1, "b"=0.1, "c"=0.1, "a/b"=0.1, "a/c"=0.1, 
                "b/c"=0.1, "a/b/c"=0.4), stateSpace)

# apply Yager's combination operator, m12 will be a massQ-object 
m12Q <- yComb(m1,m2)

# contruct a mass function from an massQ-object
m12 <- mass(m12Q)

}

