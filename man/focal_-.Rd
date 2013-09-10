\name{focal<-}
\alias{focal<-}
\title{
Replacement Functino for Focal Elements
}

\description{
Replaces focal elements of a mass function.
}
\usage{
focal(x) <- value
}

\arguments{
    \item{x}{
        a mass function
    }
    \item{value}{
        new focal elements for the mass function
    }
}

\value{
mass function with focal elements replaced.
}
\references{
Dempster, A. P. (1969), A generalization of Bayesian inference, \emph{Journal of the Royal Statistical Society}, \bold{30}, 205-247

Shafer, G., (1976), A Mathematical Theory of Evidence Princeton University Press
}

\author{
Alexander Karlsson
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass functions 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)

# replace focal elements
focal(m) <- list("a/b"=1)
}

