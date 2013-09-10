\name{mass-methods}
\docType{methods}
\alias{mass-methods}
\alias{mass,list,character-method}
\alias{mass,massQ,missing-method}

\title{Methods for Function \code{mass}}
\description{
Methods for constructing a mass function. For more detail see \code{\link{mass}}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "list", y = "character")}}{
Construct a mass functions by a named list of focal elements and a given state space
}

\item{\code{signature(x = "massQ", y = "missing")}}{
Construct a mass function from a \code{\link{massQ-class}} object
}
}}

\author{
Alexander Karlsson
}

\keyword{methods}

