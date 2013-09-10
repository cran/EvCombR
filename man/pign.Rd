\name{pign}
\alias{pign}
\title{
Pignistic Tranformation
}
\description{
The pignistic transformation transforms a mass function into a probability function. 
}
\usage{
pign(x)
}
\arguments{
  \item{x}{
    a mass function
}
}
\value{
a singleton credal set
}

\references{
Smets, P. & Kennes, R. (1994), The transferable belief model, \emph{Artificial Intelligence}, \bold{66}, 191-234
}

\author{
Alexander Karlsson
}


\seealso{
\code{\link{relPl}}
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass function 
m <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)

# obtaina singleton credal set
c <- pign(m)

}

