\name{[<--methods}
\docType{methods}
\alias{[<--methods}
\alias{[<-,credal,ANY,ANY-method}
\alias{[<-,mass,character,missing-method}

\title{Methods for Function \code{[<-}}
\description{
Replace part of an evidence structure
}
\section{Methods}{
\describe{
\item{\code{signature(x="credal", i="ANY", j="ANY", value="ANY")}}{
    Replace probabilities
}

\item{\code{signature(x="mass", i="character", j="missing", value="ANY")}}{
    Replace focal element(s)

}
}}

\author{Alexander Karlsson}

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

# construct credal sets with the given state space
c <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace) 

# replace first and second extreme point
c[1:2,] <- rbind(c(0.1, 0.1, 0.8), c(0.2, 0.2, 0.6))

# mass function 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)  

# switch mass on focal elements "b" and "a/b/c" 
temp <- m["b"]
m["b"] <- m["a/b/c"]
m["a/b/c"] <- temp

}

\keyword{methods}
