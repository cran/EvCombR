\name{[[<--methods}
\docType{methods}
\alias{[[<--methods}
\alias{[[<-,mass,character,missing-method}

\title{Methods for Function \code{[[<-}}
\description{
Replace part of an evidence structure
}
\section{Methods}{
\describe{

\item{\code{signature(x="mass", i="character", j="missing", value="ANY")}}{
    Replace focal element(s)
}
}}

\author{Alexander Karlsson}

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

# mass function 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)  

# obtain value only
m[["a"]]

}

\keyword{methods}
