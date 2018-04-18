scenario.names <- function()
{
  nms <- c("FiveYrCG21",
           "FiveYrCG703",
           "FiveYrCG710",
           "FiveYrGG22",
           "FiveYrGG704",
           "FiveYrGG7",
           "FiveYrCG22",
           "FiveYrCG704",
           "FiveYrCG7",
           "FiveYrGG701",
           "FiveYrGG707",
           "FiveYrGGPG",
           "FiveYrCG701",
           "FiveYrCG707",
           "FiveYrCGPG",
           "FiveYrGG702",
           "FiveYrGG709",
           "FiveYrCG702",
           "FiveYrCG709",
           "FiveYrGG21",
           "FiveYrGG703",
           "FiveYrGG710")
  return(nms)
}

scenario.setup <- function(nm,
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
{

  ##------ Load Dataset ---------
  dd <- paste0("data(",nm,")")
  eval(parse(text = dd))
  data <- eval(parse(text = nm))
  ##-----------------------------

  ##------- Break down data -------
  Y <- data[,c(4)]
  if(is.null(w.X))
  {
    w.X <- 5:23
  }
  X <- data[,w.X]
  if(is.null(w.Z))
  {
    w.Z <- c(25:44)
  }
  Z <- data[,w.Z]
  p <- dim(data)[2]
  if(is.null(w.W))
  {
    w.W <- c(47:52,54:p)
  }
  W <- data[,w.W]
  W <- cbind(W,1)
  p.W <- dim(W)[2]
  p.X <- dim(X)[2]
  Y.g <- cbind(X,Y)
  ##------------------------------

  ##------- Now form U -----------
  U <- list()
  for(i in 1:dim(X)[2])
  {
    U[[i]] <- cbind(W,Z)
  }
  U[[dim(X)[2] + 1]] <- cbind(W,X)
  ##-------------------------------

  ##--- Now make model parameters ---
  pi.M <- list()
  Theories <- list()
  Theories.prob <- list()
  for(j in 1:length(U))
  {
    p.U <- dim(U[[j]])[2]

    if(is.null(pi.M.fix))
    {
      pi.M[[j]] <- c(rep(1, p.W), rep(0.5, p.X))
    }else{
      pi.M[[j]] <- pi.M.fix
    }

    if(is.null(Theories.fix))
    {
      Theories[[j]] <- c(rep(1,6),rep(2,p.W - 7),3,4,5,6,6,7,8,8,8,8,9,10,11,12,12,12,13,13,13,13,14)
    }else{
      Theories[[j]] <- Theories.fix
    }
    
    m <- max(Theories[[j]])

    if(is.null(Theories.prob.fix))
    {
      Theories.prob[[j]] <- c(rep(1,3), rep(0.5, m - 3))
    }else{
      Theories.prob[[j]] <- Theories.prob.fix
    }
    if(!is.null(robustness.fix.theory))
    {
      Theories.prob[[j]][robustness.fix.theory] <- robustness.fix.prob
    }
  }
  ##-----------------------------------

  ##----- Object to return ------------
  l <- list(Y.g = Y.g,U = U, pi.M = pi.M, Theories = Theories, Theories.prob = Theories.prob,
            Z = Z, which.intercept = which.intercept)
  ##-----------------------------------

  return(l)
}

run.scenario <- function(nm, s = 1e5,
                         b = round(s/10),
                         print.every = 100,
                         out.loc = "./Results/",
                         robustness.fix.theory = NULL,
                         robustness.fix.prob = 1,
                         Theories.prob.fix = NULL,
                         Theories.fix = NULL,
                         pi.M.fix = NULL,
                         w.W = NULL,
                         w.Z = NULL,
                         which.intercept = 1,
                         w.X = NULL)
{

  print(paste("Starting scenario",nm,Sys.time()))
  
  ##--------- Get Setup ---------
  l <- scenario.setup(nm,
                      robustness.fix.theory,
                      robustness.fix.prob,
                      Theories.prob.fix = Theories.prob.fix,
                      Theories.fix = Theories.fix,
                      pi.M.fix = pi.M.fix,
                      w.W = w.W,
                      w.Z = w.Z,
                      which.intercept = which.intercept,
                      w.X = w.X)
  ##-----------------------------

  ##-------- Main Run -----------
  Results <- sur.bma.theory(l$Y.g, l$U, s=s, b=b,
                            full = FALSE,
                            pi.M = l$pi.M,
                            Theories = l$Theories,
                            Theories.prob = l$Theories.prob,
                            print.every = print.every,
                            which.intercept = l$which.intercept)
  ##------------------------------

  ##------- Write main results -----
  U.lab <- list()
  for(j in 1:length(l$U))
  {
    U.lab[[j]] <- c(colnames(l$U[[j]]))
  }
  tbl <- sur.bma.summary(Results, U.lab)
  capture.output(tbl,file=paste0(out.loc,"/Res_",nm,"_IVBMA.txt"))
  ##---------------------------------------------------

  ##--- Plot the posterior model size distribution ----
  png(paste0(out.loc,"/Res_",nm,"_PMSD.png"))
  for(j in 1:length(Results$M))
  {
    n.el <- rowSums(Results$M[[j]])
    w.n <- sort(unique(n.el))
    if(j < length(Results$M))
      {
        m <- paste("Endogenous Model",j)
      }else{
        m <- paste("Outcome Model")
      }

    plot(table(n.el)/sum(table(n.el)), type = "p",xlab="Model Size",ylab="Probability",pch = 20,main = m)
    lines(table(n.el)/sum(table(n.el)))
  }
  dev.off()
  ##---------------------------------------------

  ##---------------------------------------------
  ## Get Top models
  top.models <- list()
  for(j in 1:length(l$U))
  {
    top.models[[j]] <- model.tabulation(Results$M[[j]],U.lab[[j]],n=5)
  }
  capture.output(top.models,file=paste0(out.loc,"/Res_",nm,"_TopModels.txt"))
  ##---------------------------------------------

  ## --------- Run Diagnostics ------------------
  ## This is the probability that the instruments violate the instrument assumption
  diagn <- sur.bma.iv.diagnostics(Results,l$Z)
  p.sargan <- mean(diagn$BTIV.collapsed.cond)
  capture.output(p.sargan,file=paste0(out.loc,"/Res_", nm,"_Sargan.txt"))
  ##--------------------------------------------
  
  ##-------------------------------------------
  ## Get Theory Probabilities
  ## Return Theory probability.  
  prob.theory <- posterior.theory.probability(Results$M, l$Theories)
  capture.output(prob.theory,file=paste0(out.loc,"/Res_",nm,"_TheoriesProb.txt"))
  ##-------------------------------------------

  print(paste("Finished scenario",nm,Sys.time()))
  
  return(TRUE)
}
                 
