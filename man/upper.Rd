\name{upper}
\alias{upper}
\title{
Upper Bounds Based on Evidence Structure
}
\description{
Calculate the upper bounds for a vector of sets         
}
\usage{
upper(x, sets)
}

\arguments{
\item{x}{
credal set or mass function
}
\item{sets}{
vector of sets where each set is represented by state names separated by "/". If sets are missing, upper bounds on singletons are calculated.
}
}

\value{
upper bound of mass or probability for each set in the vector sets or if sets is missing upper bounds on singletons               
}

\references{
Shafer, G., (1976), A mathematical theory of evidence, Princeton University Press

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

# calculate upper bounds
upper(m, c("a", "a/b"))
upper(c, c("a", "a/b")) 

# upper bounds on singletons
upper(m)

}
