\name{[[-methods}
\docType{methods}
\alias{[[-methods}
\alias{[[,mass,character,missing-method}
\title{Methods for Function \code{[[}}
\description{
Methods for function \code{[[}
}
\section{Methods}{
\describe{
    \item{\code{signature(x="mass", i="character", j="missing")}}{
        Extract a single focal element from the list of focal elements
    }
}}

\author{Alexander Karlsson}

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

#mass functions 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)  

# extract focal element
m[["a"]]
}                          


\keyword{methods}

