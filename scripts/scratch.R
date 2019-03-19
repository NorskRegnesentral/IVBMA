nm = nms[1]
b = round(s/10)
print.every = 100
out.loc = "./Results/"
robustness.fix.theory = NULL
robustness.fix.prob = 1
Theories.prob.fix = NULL
Theories.fix = NULL
pi.M.fix = NULL
w.W = NULL
w.Z = NULL
which.intercept = 1
w.X = NULL

print.every = ceiling(s / 1e2)
out.loc = out.loc
Theories.prob.fix = Theories_prob
Theories.fix = Theories
pi.M.fix = pi.M
w.X = w_X
w.Z = w_Z
w.W = w_W
which.intercept = which_intercept

s=1e3
b = round(s/10)
full = FALSE,
odens = min(c(5e3,s-b))
pi.M = NULL
Theories = NULL
Theories.prob = NULL
print.every = round(s/10)
which.intercept = 1

full = FALSE
pi.M = l$pi.M
Theories = l$Theories
Theories.prob = l$Theories.prob
print.every = print.every
which.intercept = l$which.intercept
 
Y = l$Y.g





robustness.fix.theory = NULL,
robustness.fix.prob = 1,
Theories.prob.fix = NULL,
Theories.fix = NULL,
pi.M.fix.last = NULL,
Theories.prob.fix.last = NULL,
Theories.fix.last = NULL,
pi.M.fix = NULL,
w.W = NULL,
w.Z = NULL,
which.intercept = 1,
w.X = NULL)



X = read.csv("~/pkg/IVBMA/data/R8_FiveYrCG7_SSA.txt", sep="\t")
X[,"gdp"] = X[,"gdp"] / 1e12
X[,"lgdp"] = X[, "lgdp"] / 1e12
write.table(X, file =  "~/pkg/IVBMA/data/R8_FiveYrCG7_SSA.txt", row.names = FALSE, sep ="\t")

