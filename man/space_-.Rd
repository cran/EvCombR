\name{space<-}
\alias{space<-}
\title{
Replacement Function for State Space
}
\description{
Replace the names of the state space
}
\usage{
space(x) <- value
}

\arguments{
\item{x}{
    mass function or credal set
}
\item{value}{
new state space given as a character vector
}
}
\value{
new mass function or credal set with the state space replaced
}

\author{
Alexander Karlsson
}

\seealso{
\code{\link{focal<-}}
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# construct mass function 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)

# replace state space
space(m) <- c("d", "e", "f")
}
