\name{EvCombR News}
\title{News for EvCombR}
\encoding{UTF-8}

\section{Version 0.1-3 (2022-04-20)}{
    \itemize{
        \item License changed to GPL (>=3)
		\item Description and namespace files changed
		\item Added function EvCombRLicense	
    }
}      

\section{Version 0.1-2 (2014-04-22)}{
    \itemize{
        \item \code{[[<-} has been added for mass functions.  
        \item If the upper bounds of the probability intervals in \code{credal} are missing, 
              \code{1} is assumed. Note that in previous versions this 
              case was interpreted as a single probability function but now you have to construct a credal set containing only one
              probability function by specifying the same lower and upper bounds on probability intervals, e.g., 
              \code{credal(c(0.1., 0.4, 0.5), c(0.1., 0.4, 0.5), c("a", "b", "c")}.
        \item \code{lower} and \code{upper} return bounds for singleton sets if no sets are provided. 
        \item Focal sets in the mass functions are now sorted by cardinality.
        \item \code{"["}, \code{"[["}, \code{"[<-"} and \code{"[[<-"} are order invariant for mass functions 
              (i.e., \code{m["a/b"]} and \code{m["b/a"]} are equivalent). 
    }
}     

\section{Version 0.1-1 (2014-01-22)}{
    \itemize{
        \item \code{lower} and \code{upper} have been vectorized with respect to sets.  
        \item \code{lower} is now invariant to the specification of the set. 
    }
}     

\section{Version 0.1 (2013-09-10)}{
    \itemize{
        \item First release containing: Bayesian, credal (based on intervals), Yager's, Dempster's, and
              modified Dempster's combination. The package also contains function 
              for calculating pignistic and relative plausibility transforms, upper
              and lower bounds, and a discounting operator for mass functions. 
    }
}

                      
