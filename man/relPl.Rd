\name{relPl}
\alias{relPl}
\title{
Relative Plausibility Transform
}
\description{
The relative plausibility transform transform a mass function to a probability function
}
\usage{
relPl(x)
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
Cobb, B. & Shenoy, P. (2006), On the plausibility transformation for translating belief function models to probability models, \emph{International Journal of Approximate Reasoning}, \bold{42}, 3, 314 - 330
}
\author{
Alexander Karlsson
}

\seealso{
\code{\link{pign}}
}


\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass function 
m <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)

# obtaina singleton credal set
c <- relPl(m)
}

