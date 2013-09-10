\name{yComb-methods}
\docType{methods}
\alias{yComb-methods}
\alias{yComb,mass,mass-method}
\alias{yComb,list,missing-method}

\title{Methods for Function \code{yComb}}
\description{
Combine mass functions using Yager's combination operator. For more detail
see \code{\link{yComb}}.    
}
\section{Methods}{
\describe{
\item{\code{signature(x = "mass", y = "mass")}}{
    Combine two mass functions using Yager's combination operator  
}

\item{\code{signature(x = "list", y = "missing")}}{
    Combine a list of mass functions using Yager's combination operator  
}
}}


\keyword{methods}

