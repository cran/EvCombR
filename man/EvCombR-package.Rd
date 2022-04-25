\name{EvCombR-package}
\alias{EvCombR-package}
\alias{EvCombR}
\docType{package}
\title{
	EvCombR - Evidence Combination in R
}
\description{
Package for combining pieces of evidence. 
}
\details{
Implements Dempster's, Yager's, modified Dempster's, Bayesian, and credal combination (based on intervals). 
}
\author{
Alexander Karlsson

Maintainer: Alexander Karlsson <alexander.karlsson@his.se>
}
\references{
Dempster, A. P. (1969), A generalization of Bayesian inference, \emph{Journal of the Royal Statistical Society}, \bold{30}, 205-247

Shafer, G. (1976), A Mathematical Theory of Evidence Princeton University Press  

Yager, R. (1987), On the Dempster-Shafer Framework and New Combination Rules,
\emph{Information Sciences 41: 93-137}.

Fixsen, D., Mahler, R. P. S. (1997), The modified Dempster-Shafer approach to classification, \emph{IEEE Transactions on Systems, Man and Cybernetics, Part A: Systems and Humans}, \bold{27}, 96-104

Arnborg, S. (2006), Robust Bayesianism: Relation to Evidence Theory, \emph{Journal of Advances in Information Fusion}, \bold{1}, 63-74

Karlsson, A., Johansson, R., and Andler, S. F. (2011), Characterization and Empirical Evaluation of Bayesian and Credal Combination Operators, \emph{Journal of Advances in Information Fusion}, \bold{6}, 150-166

}

\keyword{ package }

\examples{
# construct a state space
stateSpace <- c("a", "b", "c")

# construct credal sets with the given state space
c1 <- credal(c(0.1, 0.1, 0.1), c(0.8, 0.8, 0.8), stateSpace)
c2 <- credal(c(0.2, 0.2, 0.2), c(0.9, 0.9, 0.9), stateSpace)

# combine the credal sets
cComb(c1, c2)

# construct mass functions
m1 <- mass(list("a"=0.1, "b"=0.1 , "c"=0.4, "a/b/c"=0.4), stateSpace)
m2 <- mass(list("a"=0.2, "b"=0.2, "c"=0.2, "a/b/c"=0.4), stateSpace)

# combine the mass functin by using Dempster's combination
dComb(m1, m2)

# Yager's combination operator
yComb(m1, m2)

# modified Dempster's combination using uniform prior
mComb(m1, m2)

}
