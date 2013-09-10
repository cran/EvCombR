\name{lower}
\alias{lower}
\title{
Lower Bounds on Evidence Structure
}
\description{
Calculate the lower bounds for a specific set of states
}
\usage{
lower(x, set)
}

\arguments{
\item{x}{
credal set or mass function
}
\item{set}{
set of states represented as character where the states is separated by "/"
}
}

\value{
lower bound of mass or probability for states in the set
}

\references{
Shafer, G., (1976), A Mathematical Theory of Evidence Princeton University Press

Walley, P. (2000), Towards a unified theory of imprecise probability, \emph{International Journal of Approximate Reasoning, \bold{24}, 125-148
}}

\author{
Alexander Karlsson
}

\note{
This is equivalent to Belief in Dempster-Shafer theory
}

\seealso{
\code{\link{upper}}
}
\examples{
# state space
stateSpace <- c("a", "b", "c")

# mass function
m <- mass(list("a"=0.1, "b"=0.1 , 
               "c"=0.4, "a/b/c"=0.4), stateSpace)

# credal set
c <- credal(c(0.1, 0.1, 0.1), 
            c(0.8, 0.8, 0.8), stateSpace)

# calculate lower bounds
lower(m, "a/b")
lower(c, "a/b")

}
