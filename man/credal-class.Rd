\name{credal-class}
\Rdversion{1.1}
\docType{class}
\alias{credal-class}


\title{Class \code{"credal"}}

\description{
Represents a credal set by a set of extreme points. For more detail see \code{\link{credal}}.
}

\section{Objects from the Class}{
Objects can be created by \code{\link{credal}}.
}

\section{Slots}{
\describe{
    \item{\code{extPoints}:}{Object of class \code{"matrix"}. Each row is an extreme point of the credal set.}
}}


\section{Methods}{
\describe{
    \item{[}{\code{signature(x = "credal")}: extract an extreme point}
    \item{[<-}{\code{signature(x = "credal")}: replace and extreme point}
    \item{cComb}{\code{signature(x = "credal", y = "credal")}: combine two credal sets }
    \item{lower}{\code{signature(x = "credal", set = "character")}: calculate the lower bound for a specific set of states }
    \item{lower}{\code{signature(x = "credal", set = "missing")}: calculate the lower bounds for all singleton states}
    \item{upper}{\code{signature(x = "credal", set = "character")}: calculate the upper bound for a specific set of states  }
    \item{upper}{\code{signature(x = "credal", set = "missing")}: calculate the upper bounds for all singleton states}
    \item{extPoints}{\code{signature(x = "credal")}: access method for the slot points}
    \item{space}{\code{signature(x = "credal")}: access method for names of singleton states}
    \item{space<-}{\code{signature(x = "credal")}: replace method for names of singleton states}
    }
}

\author{
	Alexander Karlsson
}

\keyword{classes}
