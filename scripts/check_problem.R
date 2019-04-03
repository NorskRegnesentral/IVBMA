## Specification 1
rm(list = ls())

library(ivbma)
library(parallel)

setwd("~/pkg/IVBMA/scripts/")

revision_number = 1
s = 1e2

##------ Paste Here -----------
nms <- c("R1_FiveYrGG7_1",
         "R1_FiveYrCG7_1")

nm <- nms[1]
dd <- paste0("data(",nm,")")
eval(parse(text = dd))
data <- eval(parse(text = nm))
p = dim(data)[2]

w_Y = 4
w_X <- c(5:23)
w_Z <- c(24:61)
w_W <- c(65:69,71:p)

p.W = length(w_W)
p.X = length(w_X)

pi.M <- c(rep(1, p.W), rep(0.5, p.X))
Theories <- c(rep(1,5),rep(2,p.W - 6),3,4,5,6,6,7,7,7,7,8,9,10,10,10,11,11,11,11,11,12)
m = max(Theories)
Theories_prob <- c(rep(1,3), rep(0.5, m - 3))
which_intercept = 3
##----------------------------
b = round(s/10)
out.loc = "./Results/"
robustness.fix.theory = NULL
robustness.fix.prob = 1

s = s
print.every = 1
Theories.prob.fix = Theories_prob
Theories.fix = Theories
pi.M.fix = pi.M
w.X = w_X
w.Z = w_Z
w.W = w_W
which.intercept = which_intercept
