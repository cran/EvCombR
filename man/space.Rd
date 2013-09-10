\name{space}
\alias{space}
\title{
State Space of and Evidence Structure
}
\description{
This functions returns the state space of an evidence structure.
}
\usage{
space(x)
}
\arguments{
  \item{x}{
    mass function or credal set
}
}
\value{
a character vector with the names within the state space

}

\author{
Alexander Karlsson
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# construct mass function 
m <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)

# obtain state space
space(m)

}
