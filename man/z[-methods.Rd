\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{[,credal-method}
\alias{[,mass-method}
\title{Methods for Function \code{[}}
\description{
Extract part of evidence structure \code{[} 
}

\section{Methods}{

\describe{

\item{\code{signature(x = "credal", i = "ANY", j="ANY")}}{
Extract probabilities
}

\item{\code{signature(x = "mass", i= "ANY", j="ANY")}}{
Extract focal element(s)
}
}}

\author{Alexander Karlsson}

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

# construct credal sets with the given state space
c <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace) 

# extract first and second extreme point
c[1:2,]

# mass functions 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)  

# extract focal elements
m[c("a", "a/b/c")]
}

\keyword{methods}

