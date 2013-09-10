\name{cComb-methods}
\docType{methods}
\alias{cComb,credal,credal-method}
\alias{cComb,list,missing-method}

\title{Methods for Function \code{cComb}}

\description{
Combine credal sets (based on intervals) using the credal combination operator (also known as the robust Bayesian combination operator). For more detail see \code{\link{cComb}}.     
}

\section{Methods}{
\describe{
    \item{\code{signature(x = "credal", y = "credal")}}{
        Combine two credal sets using the credal combination operator
    }
    \item{\code{signature(x = "list", y = "missing")}}{
        Combine a list of credal sets using the credal combination operator
    }
}}

\keyword{methods}

