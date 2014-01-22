# =============================================================================
# EvCombR, R Package for combining evidence, Version 0.1-1 
# Copyright (c) 2014 Alexander Karlsson
#
# License: The MIT License (MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.    
# =============================================================================   

# ============== Support ====================================================== 
# This work was supported by the Information Fusion Research Program 
# in partnership with the Swedish Knowledge Foundation under
# grant 2010-0320 (URL: http://www.infofusion.se, UMIF project)
# =============================================================================


# ============== Hooks ========================================================
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("
EvCombR - Evidence Combination in R, Version 0.1-1 
Copyright (c) 2014, Alexander Karlsson
License: MIT
")
}

# ============== Class Definitions ============================================
# credal sets (extreme points)
setClass("credal", representation(extPoints="matrix"))

# mass functions (focal elements
setClass("mass", representation(focal="list", space="character"))

# mass functoins that stores mass on the empty set
setClass("massQ", representation(qEmpty="numeric"), contains="mass")

# ============== Validity functions ===========================================
# credal validity check
setValidity("credal", 
    function(object) {
        # all extreme points should sum to one
        ifelse(all(apply(object@extPoints, 1, checkStruc)), TRUE, 
                "Not a probability function") 
    })


# mass validity check
setValidity("mass", 
    function(object) {
        # all extreme points should sum to one
        ifelse(checkStruc(as(object@focal, "numeric")), TRUE, 
                "Not a mass function")
    })

# massQ validity check
setValidity("massQ", 
    function(object) {
       # check so that qEmpty (conflict) is between zero and one
       ifelse((object@qEmpty >= 0) && (object@qEmpty <= 1), TRUE, 
                "Q-value out of range")    
    })

# ================ Generics ===================================================
# credal constructor 
setGeneric("credal", function(x, y, z) standardGeneric("credal")) 

# mass constructor
setGeneric("mass", function(x,y) standardGeneric("mass"))            

# mass, focal replacement
setGeneric("focal<-", function(x,value) standardGeneric("focal<-"))    

# state space replacement for credal set or mass function
setGeneric("space<-", function(x,value) standardGeneric("space<-"))

# obtain state space of mass function or credal set
setGeneric("space", function(x) standardGeneric("space"))

# obtain points of a credal set
setGeneric("extPoints", function(x) standardGeneric("extPoints"))

# obtain focal elements of a mass function
setGeneric("focal", function(x) standardGeneric("focal")) 

# lower bounds for a set of states 
setGeneric("lower", function(x, sets) standardGeneric("lower"))       

# upper bounds for a set of states          
setGeneric("upper", function(x, sets) standardGeneric("upper"))    


# Dempster's combination operator
setGeneric("dComb", function(x, y) standardGeneric("dComb"))

# Modified Dempster-Shafer
setGeneric("mComb", function(x, y, z) standardGeneric("mComb"))

# Credal combination operator
setGeneric("cComb", function(x, y) standardGeneric("cComb"))

# Yager's combination operator
setGeneric("yComb", function(x, y) standardGeneric("yComb"))

# pignistic transform
setGeneric("pign", function(x) standardGeneric("pign"))

# relative plausibility
setGeneric("relPl", function(x) standardGeneric("relPl"))

# discounting operators
setGeneric("disc", function(x,y) standardGeneric("disc"))

# =================== Credal Methods ==========================================
# constructor, single probability function with state names
setMethod("credal", c("numeric", "character", "missing"), 
    function(x,y) {
        extPoints <- t(as(x, "matrix"))
        colnames(extPoints) <- y
        new("credal", extPoints=extPoints)
    })


# constructor, lower/upper bounds and state names
setMethod("credal", c("numeric", "numeric", "character"), 
    function(x, y, z) { 
        extPoints <- getExtPoints(x, y)
        colnames(extPoints) <- z
        new("credal", extPoints=extPoints)
    })


          
# method for indexing probabilities 
setMethod("[", c("credal", "ANY", "ANY"), function(x,i,j) x@extPoints[i,j])



# replacement method for probability
setReplaceMethod("[", c("credal", "ANY", "ANY", "ANY"), 
                 function(x, i, j, value) {
                     x@extPoints[i,j] <- value
                     x
                 })


# replacement method for names
setReplaceMethod("space", c("credal", "character"), 
    function(x, value) {
        credal(x@extPoints, value)
    }) 

# obtain state names 
setMethod("space", "credal", function(x) {colnames(x@extPoints)})
     
# get points (in matrix form)
setMethod("extPoints", "credal", function(x) x@extPoints) 

# lower bounds, for states in the set "set"
setMethod("lower", c("credal", "character"), 
    function(x, sets) {
        f <- function(set) {
            min(apply(as.matrix(x@extPoints[,strsplit(set,"/")[[1]]]), 1, sum))
        }
        sapply(sets, f)
    })


# upper bounds, for states in the set "set" 
setMethod("upper", c("credal", "character"), 
    function(x, sets) {
        f <- function(set) {             
              max(apply(as.matrix(x@extPoints[,strsplit(set,"/")[[1]]]), 1, sum))
        }
        sapply(sets, f)  
    })  



# credal combination operator, for credal sets
setMethod("cComb", c("credal", "credal"), 
    function(x, y)  {
        # first make lists of rows
        xl <- lapply(1:nrow(x@extPoints), function(i) x@extPoints[i,])
        yl <- lapply(1:nrow(y@extPoints), function(i) y@extPoints[i,])
        
        # unnormalized combination
        xyl <- unlist(lapply(xl, function(x) lapply(yl, function(y) x*y)), 
                        recursive=FALSE)
        
        # normalize and transform to a matrix
        extPoints <- t(sapply(xyl, function(x) x/sum(x)))
        
        # check if cComb was ok
        if (!all(is.finite(extPoints))) {
            stop("cComb is not possible")    
        }
        
        # return the joint credal set 
        credal(apply(extPoints, 2, min), 
               apply(extPoints, 2, max), space(x))
    })

# credal combination operator, list of credal sets
setMethod("cComb", c("list", "missing"), 
    function(x) {
        Reduce(cComb, x)    
    })  

# ===================== Mass Methods ==========================================
# constructor, list containing set=mass and the space 
setMethod("mass", c("list", "character"), 
    function(x, y) {
        # check mass function before rounding and normalization
        if (!checkStruc(as(x, "numeric"))) {
            stop("Not a mass function")
        }
 
        # remove sets with zero mass
        res <- x[x > 0]    
    
        # construct mass function
        new("mass", focal=res, space=y)
    })
    
# constructor for massQ 
setMethod("mass", c("massQ", "missing"),
    function(x) {
        # restore mass function
        stateSpace <- paste0(x@space, collapse="/")
        if (!is.null(x@focal[[stateSpace]])) {
            x@focal[[stateSpace]] <- x@focal[[stateSpace]] - x@qEmpty
        }
        # put qEmpty on the empty set
        x@focal[["ES"]] <- x@qEmpty
        # return mass function
        mass(x@focal, x@space)
    })

# method for indexing focal elements
setMethod("[", c("mass", "ANY"), function(x,i) x@focal[i])

# method for obtaining single focal element values 
setMethod("[[", c("mass", "ANY"), function(x,i) x@focal[[i]])

# obtain state names
setMethod("space", "mass", function(x) x@space)

# obtain focal sets
setMethod("focal", "mass", function(x) x@focal)         

# replacement method for focal elements
setReplaceMethod("[", c("mass", "ANY", "ANY"), 
    function(x,i,value) {
        x@focal[i] <- value
        x
    })   

                      
# replacement function for state names
setReplaceMethod("space", c("mass", "character"), 
    function(x, value) {
        mass(x@focal, value)
    })

# replacement function for focal elements
setReplaceMethod("focal", c("mass", "list"), 
                 function(x, value) {
                     mass(value, x@space)
                 })  


# lower bound (also known as Belief)
setMethod("lower", c("mass", "character"),
    function(x, sets) {
        ff <- function(set) {
            # make sets invariant to order
            f <- function(set) {
                paste0(sort(strsplit(set, "/")[[1]]), collapse="/")            
            }
            setInvar <- f(set)
            names(x@focal) <- sapply(names(x@focal), f)                           
              
            # sum over all subsets of "set" 
            sum(unlist(x@focal[powerSet(strsplit(setInvar, "/")[[1]])]),
                na.rm=TRUE)
        }          
        sapply(sets, ff)
    }) 



# upper bound (also known as plausibility)
setMethod("upper", c("mass", "character"),
    function(x, sets) {
        ff <- function(set) {        
            # returns true whenever a state within "states" is found in "mStates"
            f <- function(mStates, states) !all(is.na(match(states, mStates)))
            # apply over each focal element (set)
            i <- unlist(lapply(strsplit(names(x@focal), "/"), f, 
                               strsplit(set, "/")[[1]]))
            # return plausibility
            sum(unlist(x@focal)[i])
        }
        sapply(sets, ff)
    })  


# Dempster's combination operator, two mass functions
setMethod("dComb", c("mass", "mass"), 
    function(x, y) {
        arithOp <- function(i, j) {
            x@focal[[i]] * y@focal[[j]]    
        }
        mass(comb(x, y, intersect, Vectorize(arithOp), TRUE), x@space)
    })

# Dempster's combination operator, list of mass functions
setMethod("dComb", c("list", "missing"), 
    function(x) {
        Reduce(dComb, x)    
    })  

# Modified Dempster-Shafer, two mass functions and q
setMethod("mComb", c("mass", "mass", "list"), 
    function(x, y, z) {
        # help function for calculating the probability for a set of states 
        # using the the list z (a probability distribution) 
        q <- function(cVec) {sum(unlist(z[cVec]))}  
        
        # z is a distribution over state names
        # define the arithmetic operation
        arithOp <- function(i, j) {
            op1StateNames <- strsplit(names(x@focal), "/")[[i]]
            op2StateNames <- strsplit(names(y@focal), "/")[[j]]
            intersectStateNames <- intersect(op1StateNames, op2StateNames)
            x@focal[[i]] * y@focal[[j]] * 
                (q(intersectStateNames) / (q(op1StateNames) * q(op2StateNames)))
        }       
        # return joint evidence
        mass(comb(x, y, intersect, Vectorize(arithOp), TRUE), x@space)
    })

# Modified Dempster-Shafer, two mass functions with default uniform q
setMethod("mComb", c("mass", "mass", "missing"), 
    function(x, y) {
        # construct uniform q-function
        l <- as((rep(1/length(x@space), length(x@space))), "list")
        names(l) <- x@space        
        
        # combine
        mComb(x, y, l)
    })


# Modified Dempster-Shafer, list of mass functions and q
setMethod("mComb", c("list", "list", "missing"), 
    # contruct a two argument wrapper for the higher order function "Reduce"
    function(x, y) {
        mCombWrapper <- function(op1, op2) {
            mComb(op1, op2, y)            
        }    
        Reduce(mCombWrapper, x)    
    })

# Modified Dempster-Shafer, list of mass functions with default q (uniform)
setMethod("mComb", c("list", "missing", "missing"), 
    # contruct a two argument wrapper for the higher order function "Reduce"
    function(x) {
        Reduce(mComb, x)    
    })

# Yager's combination operator (quasi-associative), two mass functions and
# state space
setMethod("yComb", c("mass", "mass"), 
    function(x,y) {
        # check if x or y are of type massQ and if so transform to a mass 
        # function
        if (is(x, "massQ")) {
            x <- mass(x)    
        }
        if (is(y, "massQ")) {
            y <- mass(y)    
        } 
        
        # arithmetic operation for combination
        arithOp <- function(i, j) {
            x@focal[[i]] * y@focal[[j]]    
        }
        
        # unnormalized combination
        focal <- comb(x, y, intersect, Vectorize(arithOp), FALSE) 
        
        # check if there is any mass on the empty set
        if (!is(focal[["ES"]], "NULL")) {
            # obtain mass on empty set
            qEmpty <- focal[["ES"]]
            # add this mass to ignorance
            stateSpace <- paste0(x@space, collapse="/")
            focal[[stateSpace]] <- ifelse(is(focal[[stateSpace]], "NULL"), 0, 
                                            focal[[stateSpace]]) + qEmpty
            # remove empty set
            focal[["ES"]] <- NULL
        } else {
            qEmpty <- 0    
        }
        
        # construct a mass function with mass on empty set added to ignorance 
        # and that keep track of the amount of mass on the empty set
        new("massQ", mass(focal, x@space), qEmpty=qEmpty) 
    })

# Yager's combination operator (quasi-associative), list of mass functions and 
# state space
setMethod("yComb", c("list", "missing"), 
    function(x) {
        Reduce(yComb, x)       
    })

# Pignistic transformation for mass functions
setMethod("pign", "mass", 
    function(x) {
        # obtain singletons
        elem <- factor(unlist(strsplit(names(x@focal), "/")))
        # distribute mass over the singletons 
        f <- function(x, y) {
            nElements <- length(unlist(strsplit(y, "/")))
            rep(x/nElements, nElements)
        }
        # gather all contributions (from all focal sets)
        mElem <- unlist(mapply(f, x@focal, names(x@focal)))
        # return a singleton credal set
        temp <- tapply(mElem, elem, sum)
        res <- rep(0, length(x@space))
        res[match(names(temp), x@space)] <- temp
        credal(as(res, "vector"), x@space)       
})           

# Relative plausibility transformation for mass functions
setMethod("relPl", "mass",            
    function(x) {
        # obtain plausabilities
        pl <- sapply(x@space, function(set) upper(x,set))
        # obtain normalization constant
        K <- sum(pl)
        # return relative plausibility for the specified set
        credal((pl/K), x@space)
    })
       

# discounting operator for mass functions
setMethod("disc", c("mass", "numeric"),
    function(x, y) {
        # perform discounting
        x@focal <- lapply(x@focal, function(z) y*z)
        # fix discounting on state space 
        stateSpace <- paste0(x@space, collapse="/")
        x@focal[[stateSpace]] <- 1 - y + ifelse(is.null(x@focal[[stateSpace]]), 
                                            0, x[[stateSpace]])
        # return discounted mass function
        x
    })



# ========== Help functions ===================================================
# Returns the powerset    
powerSet <- function(stateNames) {
    # construct powerset
    temp <- lapply(1:length(stateNames), function(n) combn(stateNames, n))
    # return the powerset
    unlist(lapply(temp, function(x) 
                       apply(x, 2, function(x) paste0(x, collapse = "/"))))
}
      
# Obtain extreme points based on lower and upper bounds
# see paper: Probability Intervals: A Tool for Uncertain Reasoning
getExtPoints <- function(low, up) {
    # To be a matrix of extreme points
    extMat <- NULL
    
    # help function used in the recursion
    g <- function(i, probVec, lambda, indices) {
        if (lambda <= up[i] - low[i]) {
            v <- probVec[i]
            probVec[i] <- probVec[i] + lambda
            # add extreme point (if it does not already exists)
            if (!is.null(extMat)) { 
                if (!any(apply(extMat, 1, 
                	function(x) isTRUE(all.equal(x, probVec))))) {
                	extMat <<- rbind(extMat, probVec, deparse.level=0)
                }
            } else {
               # first extreme point found
               extMat <<- rbind(extMat, probVec, deparse.level=0)  
            }
            probVec[i] <- v
        } else {
            v <- probVec[i]
            probVec[i] <- up[i]
            f(probVec, lambda - up[i] + low[i], indices[-i])
            probVec[i] <- v
        }
    }   

    # recursive function to obtain the extreme points
    f <- function(probVec, lambda, indices) {
        if (length(indices) != 0) {
            sapply(indices, g, probVec, lambda, indices) 
        }
    }
    
    # start the recursion
    f(low, 1 - sum(low), 1:length(low)) 
    
    # return the extreme points
    extMat
}

# higher order function for combining evidences in the form of mass functions
comb <- function(m1, m2, setOp, arithOp, norm) {
    # get intersections
    temp <- outer(strsplit(names(m1@focal), "/"),
                  strsplit(names(m2@focal), "/"), 
                  Vectorize(function(x,y) { 
                      res <- setOp(x,y) 
                      if (length(res) == 0) {
                          "ES"
                      } else {
                          res
                      }
                  }, SIMPLIFY=FALSE))
    # tranform to factors (makes calculations easier)
    sets <- factor(sapply(temp, function(x) paste0(x, collapse = "/")))
    # get the product over the intersections
    m <- as(outer(1:length(m1@focal), 1:length(m2@focal), arithOp), "vector")
    # sum over the sets
    res <- as.list(tapply(m, sets, sum))
    # normalize 
    if (norm) {
        # remove the emptyset (if present)
        res[["ES"]] <- NULL
        # sum over the remaining focal elements
        normConst <- sum(as(res, "numeric"))
        # normalize
        res <- lapply(res, function(x) x / normConst) 
    }
    # return result (sorted by the cardinality of the sets)
    order <- sapply(strsplit(names(res), "/"), Vectorize(length))
    res[unlist(sapply(1:max(order), function(x, order) 
            {which(x==order)}, order))]
}                          



# Validity checking of evidence structures
checkStruc <- function(x) {
    # sum to one and >= 0
    isTRUE(all.equal(1, sum(x))) && all(x >= 0)
}


