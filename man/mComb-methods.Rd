\name{mComb-methods}
\docType{methods}
\alias{mComb-methods}
\alias{mComb,list,list,missing-method}
\alias{mComb,list,missing,missing-method}
\alias{mComb,mass,mass,list-method}
\alias{mComb,mass,mass,missing-method}
\title{Methods for Function \code{mComb}}

\description{
Combine mass functions using modified Dempster's combination operator. For more detail
see \code{\link{mComb}}.
}

\section{Methods}{
\describe{
\item{\code{signature(x = "mass", y = "mass", z = "list")}}{
    Combine two mass functions using modified Dempster's combination operator and a prior    
} 
\item{\code{signature(x = "mass", y = "mass", z = "missing")}}{
    Combine two mass functions using modified Dempster's combination operator and the uniform prior
} 
\item{\code{signature(x = "list", y = "list", z = "missing")}}{
    Combine a list of mass functions using modified Dempster's combination operator and a prior      
}
\item{\code{signature(x = "list", y = "missing", z = "missing")}}{
    Combine a list of mass functions using modified Dempster's combination operator and the uniform prior
}
}}

\keyword{methods}

