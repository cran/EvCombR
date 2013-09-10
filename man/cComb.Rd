\name{cComb}
\alias{cComb}


\title{
Credal Combination Operator (restricted to intervals) 
}

\description{
Combine evidence in the form of credal sets (based on intervals) using the credal combination operator (also known as the robust Bayesian combination operator). The resulting credal set is approximated by using probability intervals.
}

\usage{
cComb(x,y)
}   

\arguments{
  \item{x}{ 
  	  credal set or a list of credal sets  
}
  \item{y}{ 
  	  credal set if \code{x} is a credal set, otherwise missing
}
}

\value{
credal set   

}
\references{

Levi, I. (1983), The enterprise of knowledge, The MIT press

Arnborg, S. (2006), Robust Bayesianism: Relation to Evidence Theory, \emph{Journal of Advances in Information Fusion}, \bold{1}, 63-74

Karlsson, A., Johansson, R., and Andler, S. F. (2011), Characterization and Empirical Evaluation of Bayesian and Credal Combination Operators, \emph{Journal of Advances in Information Fusion}, \bold{6}, 150-166

}
\author{
Alexander Karlsson
}


\seealso{
\code{\link{dComb}, \link{yComb}, \link{mComb}}
}

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

# construct credal sets with the given state space
c1 <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace)
c2 <- credal(c(0.2, 0.2, 0.2), c(0.9, 0.9, 0.9), stateSpace)

# combine the credal sets
cComb(c1, c2)
# or by
cComb(list(c1, c2))

}

