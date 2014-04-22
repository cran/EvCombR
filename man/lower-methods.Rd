\name{lower-methods}
\docType{methods}
\alias{lower-methods}
\alias{lower,credal,missing-method}
\alias{lower,mass,missing-method}
\alias{lower,credal,character-method}
\alias{lower,mass,character-method}

\title{Methods for Function \code{lower}}
\description{
Calculate lower bounds for a vector of sets with respect to the evidence structure. For more detail see \code{\link{lower}}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "credal", sets = "character")}}{
obtain lower bounds for a vector of sets
}

\item{\code{signature(x = "credal", sets = "missing")}}{
obtain lower bounds for all singleton states
}

\item{\code{signature(x = "mass", sets = "character")}}{
obtain the belief, or lower bounds, for a vector of sets
}

\item{\code{signature(x = "mass", sets = "missing")}}{
obtain the belief, or lower bounds, for all singleton states
}

}}

\keyword{methods}

