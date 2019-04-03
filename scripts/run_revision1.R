## Specification 1
rm(list = ls())

library(ivbma)
library(parallel)

setwd("~/pkg/IVBMA/scripts/")

revision_number = 1
s = 1e5

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

p.W = length(w_W) + 1
p.X = length(w_X)
p.Z = length(w_Z)

pi.M = list()
for(j in 1:length(w_X)){
    pi.M[[j]] <- c(rep(1, p.W), rep(0.5, p.Z))
}
pi.M[[length(pi.M) + 1]] = c(rep(1,p.W), rep(0.5, p.X))

Theories = list()
for(j in 1:length(w_X)){
    Theories[[j]] = c(rep(1,5),rep(2,p.W - 6),3,4,5,6,6,7,7,7,7,8,9,10,10,10,11,11,11,11,11,12,4,5,6,6,7,7,7,7,8,9,10,10,10,11,11,11,11,11,12)

}
Theories[[length(w_X) + 1]] =  c(rep(1,5),rep(2,p.W - 6),3,4,5,6,6,7,7,7,7,8,9,10,10,10,11,11,11,11,11,12)

Theories_prob = list()

for(j in 1:length(Theories)){
    m = max(Theories[[j]])
    Theories_prob[[j]] <- c(rep(1,3), rep(0.5, m - 3))
}

which_intercept = 3
##----------------------------

helper <- function(j)
{
  out.loc <- paste0("./Results/Revision_",revision_number, "_Run_",j,"/")
  if(!file.exists(out.loc)) dir.create(out.loc)
  run.scenario(nms[j],
               s = s,
               print.every = ceiling(s / 1e2),
               out.loc = out.loc,
               Theories.prob.fix = Theories_prob,
               Theories.fix = Theories,
               pi.M.fix = pi.M,
               w.X = w_X,
               w.Z = w_Z,
               w.W = w_W,
               which.intercept = which_intercept)
}

l = mclapply(1:2, "helper", mc.cores = 2, mc.silent = FALSE)
    
quit(save = "no")
               

