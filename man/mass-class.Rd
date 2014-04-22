\name{mass-class}
\Rdversion{1.1}
\docType{class}
\alias{mass-class}


\title{Class \code{"mass"}}

\description{
Represents a mass function by a list of focal elements and corresponding mass. For more detail see \code{\link{mass}}.
}

\section{Objects from the Class}{
Objects can be created by \code{\link{credal}}.
}

\section{Slots}{
\describe{
    \item{\code{focal}:}{a list of focal elements represented by statenames seperated by "/"}
    \item{\code{space}:}{the state space represented by a character vector}
}}

\section{Methods}{
\describe{
    \item{[}{\code{signature(x = "mass", i = "character", j = "missing")}: extract focal elements}
    \item{[[}{\code{signature(x = "mass", i = "character", j = "missing")}: extract a single focal element}
    \item{[<-}{\code{signature(x="mass", i="character", j="missing", value="ANY")}: replace focal elements}    
    \item{[[<-}{\code{signature(x="mass", i="character", j="missing", value="ANY")}: replace a single focal element}
    \item{dComb}{\code{signature(x = "mass", y = "mass")}: combine two mass functions by Dempster's combination}
    \item{focal}{\code{signature(x = "mass")}: access focal elements}
    \item{focal<-}{\code{signature(x = "mass")}: replace focal elements }
    \item{lower}{\code{signature(x = "mass", set = "character")}: calculate the lower bounds for some focal element}
     \item{lower}{\code{signature(x = "mass", set = "missing")}: calculate the lower bounds for singletons}   
    \item{mComb}{\code{signature(x = "mass", y = "mass", z = "function")}: combine two mass functions by modified Dempster's combination using a prior distribution z }
    \item{mComb}{\code{signature(x = "mass", y = "mass", z = "missing")}: combine two mass functions by modified Dempster's combination using a uniform prior distribution z}
    \item{pign}{\code{signature(x = "mass")}: calculate the pignistic transformation for single states}
    \item{relPl}{\code{signature(x = "mass")}: calculate the relative plausibility for single states}
    \item{space}{\code{signature(x = "mass")}: access the state space (frame of discernment)}
    \item{space<-}{\code{signature(x = "mass")}: replace the state space (frame of discernment)}
    \item{upper}{\code{signature(x = "mass", set = "character")}: calculate the upper bound for some focal element }
    \item{upper}{\code{signature(x = "mass", set = "character")}: calculate the upper bounds for singletons}    
    \item{yComb}{\code{signature(x = "mass", y = "mass")}: combine two mass functions using Yager's rule}
	\item{disc}{\code{signature(x = "mass", y = "numeric")}: discount mass function}
	}
}

\references{
Dempster, A. P. (1969), A generalization of Bayesian inference, \emph{Journal of the Royal Statistical Society}, \bold{30}, 205-247

Shafer, G., (1976), A Mathematical Theory of Evidence Princeton University Press

Yager, R. (1987), On the Dempster-Shafer Framework and New Combination Rules,
\emph{Information Sciences 41: 93-137}.     

Fixsen, D., Mahler, R. P. S. (1997), The modified Dempster-Shafer approach to classification, \emph{IEEE Transactions on Systems, Man and Cybernetics, Part A: Systems and Humans}, \bold{27}, 96-104   

}

\author{
Alexander Karlsson
}


\keyword{classes}
