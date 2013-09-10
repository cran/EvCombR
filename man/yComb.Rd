\name{yComb}
\alias{yComb}
\title{
Yager's Combination Operator
}
\description{
Combine evidence in the form of mass functions using Yager's combination operator. 
}
\usage{
yComb(x,y)
}

\arguments{
\item{x}{
    single mass function or a list of mass functions      
}
  \item{y}{
    single mass function if \code{x} is a single mass function, otherwise missing     
}
}

\value{
mass function (\code{\link{massQ-class}}) 
}

\note{
Yager's combination operator is quasi-associative and therefore we need to keep track of the mass on the empty set by using the class \code{massQ}.   

}

\references{
Yager, R. (1987), On the Dempster-Shafer Framework and New Combination Rules,
\emph{Information Sciences 41: 93-137}.    
}

\author{
Alexander Karlsson
}


\seealso{
\code{\link{dComb}, \link{mComb}, \link{cComb}}
} 

\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass functions 
m1 <- mass(list("a"=0.1, "a/b/c"=0.9), stateSpace)
m2 <- mass(list("b"=0.2, "a/b/c"=0.8), stateSpace)

# Yager's combination
yComb(m1, m2)
# or 
yComb(list(m1, m2)) 

}



