\name{upper-methods}
\docType{methods}
\alias{upper-methods}
\alias{upper,credal,character-method}
\alias{upper,credal,missing-method}
\alias{upper,mass,character-method}
\title{Methods for Function \code{upper}}
\description{
Calculate lower bounds for a vector of sets with respect to the evidence structure. For more detail see \code{\link{upper}}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "credal", sets = "character")}}{
obtain upper bounds for a vector of sets 
}


\item{\code{signature(x = "mass", sets = "character")}}{
obtain the plausability, or upper bounds, for a vector of sets     
}
}}
\keyword{methods}
