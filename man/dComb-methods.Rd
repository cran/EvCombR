\name{dComb-methods}
\docType{methods}
\alias{dComb-methods}
\alias{dComb,mass,mass-method}
\alias{dComb,list,missing-method}

\title{Methods for Function \code{dComb}}

\description{
Combine mass functions using Dempster's combination operator. For more detail
see \code{\link{dComb}}.
}

\section{Methods}{
\describe{
\item{\code{signature(x = "mass", y = "mass")}}{
    Combine two mass functions using Dempster's combination operator
} 
\item{\code{signature(x = "list", y = "missing")}}{
    Combine a list of mass functions using Dempster's combination operator
}
}}

\author{Alexander Karlsson}

\keyword{methods}
