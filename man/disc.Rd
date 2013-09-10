\name{disc}
\alias{disc}
\title{
Discounting Operator
}
\description{
Discounts a mass function.
}
\usage{
disc(x,y)
}

\arguments{
\item{x}{
a mass function
}
\item{y}{
degree of reliability
}
}
\value{
mass function
}

\references{
Smets, P. (2000), Data Fusion in the Transferable Belief Model, Proceedings of the Third International Conference on Information Fusion  
}

\author{
Alexander Karlsson
}


\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass function 
m <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)

# source is only 80% reliable
mDisc <- disc(m, 0.8)
}

