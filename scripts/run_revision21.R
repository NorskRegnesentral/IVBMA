## Specification 21
rm(list = ls())

library(ivbma)
library(parallel)

setwd("~/pkg/IVBMA/scripts/")

revision_number = 21
s <- 1e2

##------ Paste Here -----------
nms <- c("R10_FiveYrCG7_NN8")
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
    w_W <- c(51:53,45:50,54:p)

    p.W = length(w_W) + 1
    p.X = length(w_X)
    
    pi.M <- c(rep(0.5,3), rep(1, p.W - 3), rep(0.5, p.X))
    Theories <- c(rep(1,3),rep(2,6),rep(3,p.W - 10),4,5,6,7,7,8,8,8,8,9,10,11,11,11,12,12,12,12,12,13)
    m = max(Theories)
    Theories_prob <- c(0.5, rep(1,3), rep(0.5, m - 4))
    which_intercept = 4

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
