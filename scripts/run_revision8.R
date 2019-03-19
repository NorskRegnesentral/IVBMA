## Specification 8
rm(list = ls())

library(ivbma)
library(parallel)

setwd("~/pkg/IVBMA/scripts/")

revision_number = 8
s <- 1e5

##------ Paste Here -----------
nms <- c("R7_FiveYrGG7_1",
           "R7_FiveYrCG7_1",
           "R7_FiveYrGG7_2",
           "R7_FiveYrCG7_2")
##----------------------------


helper <- function(j)
{
    nm <- nms[j]
    dd <- paste0("data(",nm,")")
    eval(parse(text = dd))
    data <- eval(parse(text = nm))
    p = dim(data)[2]

    w_Y = 4
    w_X <- c(5:23)
    w_Z <- c(24:42)
    w_W <- c(45:50,52:p)

    p.W = length(w_W)
    p.X = length(w_X)
    
    pi.M <- c(rep(1, p.W), rep(0.5, p.X))
    Theories <- c(rep(1,6),rep(2,p.W - 7),3,4,5,6,6,7,7,7,7,8,9,10,10,10,11,11,11,11,11,12)
    m = max(Theories)
    Theories_prob <- c(rep(1,3), rep(0.5, m - 3))
    which_intercept = 3

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

l = mclapply(1:length(nms), "helper", mc.cores = length(nms), mc.silent = FALSE)

quit(save = "no")
