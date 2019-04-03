## Specification 7
rm(list = ls())

library(ivbma)
library(parallel)

setwd("~/pkg/IVBMA/scripts/")

revision_number = 7
s <- 1e5

##------ Paste Here -----------
nms <- c("R6_FiveYrCG7_DebtHigh",
         "R6_FiveYrCG7_DebtLow",
         "R6_FiveYrCG7_DemHigh",
         "R6_FiveYrCG7_DemLow",
         "R6_FiveYrCG7_gdpHigh",
         "R6_FiveYrCG7_gdpLow",
         "R6_FiveYrCG7_IneqHigh",
         "R6_FiveYrCG7_IneqLow",
         "R6_FiveYrCG7_PolCHigh",
         "R6_FiveYrCG7_PolCLow",
         "R6_FiveYrCG7_PolRHigh",
         "R6_FiveYrCG7_PolRLow",
         "R6_FiveYrCG7_PopHigh",
         "R6_FiveYrCG7_PopLow",
         "R6_FiveYrCG7_TradeHigh",
         "R6_FiveYrCG7_TradeLow",
         "R6_FiveYrGG7_DebtHigh",
         "R6_FiveYrGG7_DebtLow",
         "R6_FiveYrGG7_DemHigh",
         "R6_FiveYrGG7_DemLow",
         "R6_FiveYrGG7_gdpHigh",
         "R6_FiveYrGG7_gdpLow",
         "R6_FiveYrGG7_IneqHigh",
         "R6_FiveYrGG7_IneqLow",
         "R6_FiveYrGG7_PolCHigh",
         "R6_FiveYrGG7_PolCLow",
         "R6_FiveYrGG7_PolRHigh",
         "R6_FiveYrGG7_PolRLow",
         "R6_FiveYrGG7_PopHigh",
         "R6_FiveYrGG7_PopLow",
         "R6_FiveYrGG7_TradeHigh",
         "R6_FiveYrGG7_TradeLow")


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

    p.W = length(w_W) + 1
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

