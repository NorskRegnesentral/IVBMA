## This short script shows how to set up the routines in a simple use case
rm(list = ls())
 
library(ivbma)
library(parallel)

## Load data.  You might replace this with a read.csv, etc., as necessary
data(FiveYrGG7)
data = FiveYrGG7

## Select which columns are response, endogenous, instrument or exogenous variables
p <- dim(data)[2]
w.Y = 4 ## What column is the Y (response) variable
w.X <- 5:23 ## What are the X (endogenous) variables
w.Z <- c(24:42) ## What are the Z (instrumental) variables)
w.W <- c(45:50, 52:p) ## What are the W (exogenous) variables

## Break down data
Y <- data[,c(4)]
X <- data[,w.X]
Z <- data[,w.Z]
W <- data[,w.W]
W <- cbind(W,1) ## remember to add an intercept to the W matrix
Y.g <- cbind(X,Y) ## The full matrix of LHS variables
p.W <- dim(W)[2]
p.X <- dim(X)[2]

## Form the collection of RHS variables
U <- list()
for(i in 1:dim(X)[2])
{
    U[[i]] <- cbind(W,Z)
}
U[[dim(X)[2] + 1]] <- cbind(W,X)

## Set up the prior proabability for inclusion of each variable
pi.M <- list()
for(j in 1:length(U)){
    pi.M[[j]] = rep(0.5, dim(U[[j]])[2])
}

## Set up the theories and theory level probabilities
## Right now each variable is it's own theory and that theory
## has priori inclusion probability 0.5
Theories <- list()
Theories.prob <- list()
for(j in 1:length(U)){
    Theories[[j]] = 1:dim(U[[j]])[2]
    Theories.prob[[j]] = rep(0.5, dim(U[[j]])[2])
}


##--- Now run main routine ---
s = 1e4 ## Iterations
b = 1e3 ## Burn in
Results = sur.bma.theory(Y.g, U, s, b,
                         pi.M = pi.M,
                         Theories = Theories,
                         Theories.prob = Theories.prob,
                         which.intercept = dim(W)[2], ## This makes sure the intercept gets included
                         print.every = 1)
                         

                         
