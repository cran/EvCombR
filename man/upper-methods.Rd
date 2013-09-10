\name{upper-methods}
\docType{methods}
\alias{upper-methods}
\alias{upper,credal,character-method}
\alias{upper,credal,missing-method}
\alias{upper,mass,character-method}
\title{Methods for Function \code{upper}}
\description{
Calculate upper bounds for some states with respect to the evidence structure. For more detail see \code{\link{upper}}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "credal", set = "character")}}{
 obtain the lower bound for a specific set of states   
}


\item{\code{signature(x = "mass", set = "character")}}{
obtain the Belief, or lower bounds, for a specific set of states
}
}}
\keyword{methods}
