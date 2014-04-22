\name{credal-methods}
\docType{methods}
\alias{credal-methods}
\alias{credal,matrix,character,missing-method}
\alias{credal,numeric,missing,character-method}
\alias{credal,numeric,numeric,character-method}
\title{Methods for Function \code{credal}}
\description{
Methods for constructing a credal set. For more detail see \code{\link{credal}}. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "numeric", y = "missing", z = "character")}}{
Construct a credal set based on the lower bounds of probability intervals for states (\code{1} will be the upper bound for all probability intervals)
}

\item{\code{signature(x = "numeric", y = "numeric", z = "character")}}{
Construct a credal based on probability intervals for states
}
}}

\author{
Alexander Karlsson
}

\keyword{methods}




