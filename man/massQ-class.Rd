\name{massQ-class}
\Rdversion{1.1}
\docType{class}
\alias{massQ-class}


\title{Class \code{"massQ"}}
\description{
Class that maintains information about the mass on the empty set. The class is used for Yager's combination operator 
}
\section{Objects from the Class}{
A \code{massQ}-object is obtained as a result of Yager's combination operator \code{\link{yComb}}. 

}
\section{Slots}{
  \describe{
    \item{\code{qEmpty}:}{mass on the empty set with respect to the previous combination}
    \item{\code{focal}:}{a list of focal elements represented by statenames seperated by "/"}
    \item{\code{space}:}{the state space represented by a character vector}
  }
}
\section{Extends}{
Class \code{"\linkS4class{mass}"}, directly.
}
\section{Methods}{
  All methods inherited from \code{\link{mass-class}} and in addition:
  \describe{
    \item{mass}{\code{signature(x = "massQ", y = "missing")}: convert the \code{massQ}-object to a \code{mass}-object}
	 }
}
\references{
Yager, R. (1987), On the Dempster-Shafer Framework and New Combination Rules,
\emph{Information Sciences 41: 93-137}.     
}

\author{
Alexander Karlsson
}

