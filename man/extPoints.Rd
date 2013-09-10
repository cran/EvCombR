\name{extPoints}
\alias{extPoints}
\title{
Extreme Points of a Credal Set
}
\description{
Returns the extreme points of a credal set 
}
\usage{
extPoints(x)
}
\arguments{
  \item{x}{
    a credal set
}
}
\value{
a matrix where the extreme points are stored by row
}


\author{
Alexander Karlsson
}


\seealso{
\code{\link{lower}}, \code{\link{upper}}
}

\examples{
# state space
stateSpace <- c("a", "b", "c")

# construct credal set
c <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace)

# obtain extrem points
eMat <- extPoints(c)

}

