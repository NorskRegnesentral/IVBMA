## Gets the fitted value for Y for a given equation r, possibly excluding the theory t.
sur.bma.theory.residual <- function(theta,D,r,exclude.t = NULL)
{
  Y <- D$Y[,r]
  n.t <- D$n.t[r]
  K <- theta$K
  for(s in 1:D$R)
  {
    if(s != r)
    {
      Y <- Y + K[s,r]/K[r,r] * (D$Y[,s] - D$U[[s]] %*% theta$beta[[s]])
    }
  }

  if(is.null(exclude.t))return(Y)

  w.exclude <- D$w.t[[r]][[exclude.t]]
  beta.temp <- theta$beta[[r]]
  beta.temp[w.exclude] <- 0
  Y <- Y - D$U[[r]] %*% beta.temp
  return(Y)
}

sur.bma.theory.update <- function(theta,D)
{

  n <- D$n
  K <- theta$K
  ##-------- Get residuals --------
  eps <- matrix(0, D$n, D$R)
  ##-------------------------------

  ##----------- Loop through equations ---------
  for(r in 1:D$R)
  {

    w.r <- D$w.t[[r]][[ D$which.intercept[r] ]] ## Everything in the intercept theory is in

    ##-------- Begin looping through all theories -----
    if(D$n.t[r] > 1)
      {
        for(t in 1:D$n.t[r])
        {

          if( !(t == D$which.intercept[r]) )
          {
            Y.residual <- sur.bma.theory.residual(theta,D,r,exclude.t=t)
            w.t <- D$w.t[[r]][[t]]
            if(!all(D$pi.M[[r]][w.t] == 1))
            {
              ##-------- Extract ------------
              M.curr <- theta$M[[r]][[t]]
              U.t <- D$U[[r]][, w.t, drop=FALSE]
              p.U <- dim(U.t)[2]
              ##-----------------------------
              
              ##------ Update Model ---------
              pp <- (1 - M.curr) * D$pi.M[[r]][w.t] + M.curr * (1 - D$pi.M[[r]][w.t])
              w <- sample(1:p.U,1,prob= pp)
              
              M.new <- M.curr
              M.new[w] <- 1 - M.new[w]
              if(any(M.new == 1))
              {
                
                ##------- Relevant calculations for new model -----
                w.1 <- which(M.curr == 1)
                p.1 <- sum(M.curr)
                U.1 <- U.t[, w.1, drop = FALSE]
                Omega.1 <- diag(p.1) + K[r, r] * t(U.1) %*% U.1
                lambda.1 <- K[r, r] * t(Y.residual) %*% U.1 %*% solve(Omega.1)
                post.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1) ## Integrated probability
                prior.1 <- mydet(D$C[[r]][[t]][which(M.curr == 1), which(M.curr == 1)])
                score.1 <- post.1 + prior.1
                ##-------------------------------------------------
                
                ##------- Relevant calculations for new model -----
                w.2 <- which(M.new == 1)
                p.2 <- sum(M.new)
                U.2 <- U.t[, w.2, drop = FALSE]
                Omega.2 <- diag(p.2) + K[r, r] * t(U.2) %*% U.2
                lambda.2 <- K[r, r] * t(Y.residual) %*% U.2 %*% solve(Omega.2)
                post.2 <- 0.5 * lambda.2 %*% Omega.2 %*% t(lambda.2) - 0.5 * mydet(Omega.2)
                prior.2 <- mydet(D$C[[r]][[t]][which(M.new == 1), which(M.new == 1)])
                score.2 <- post.2 + prior.2
                ##-------------------------------------------------
                
                ##------- Make choice --------------
                alpha <- score.2 - score.1
                if (log(runif(1)) < alpha)
                {
                  M.curr <- M.new
                  Omega.1 <- Omega.2
                  lambda.1 <- lambda.2
                  post.1 <- post.2
                }
                ##----------------------------------
              }else{
                w.1 <- which(M.curr == 1)
                p.1 <- sum(M.curr)
                U.1 <- U.t[, w.1, drop = FALSE]
                Omega.1 <- diag(p.1) + K[r, r] * t(U.1) %*% U.1
                lambda.1 <- K[r, r] * t(Y.residual) %*% U.1 %*% solve(Omega.1)
                post.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1) ## Integrated probability
              }        
              ##------- Finished updating model ----
            }else{
              ##----- If this is a model with everything included in the prior --------
              M.curr <- theta$M[[r]][[t]]
              w.1 <- which(M.curr == 1)
              p.1 <- sum(M.curr)
              U.t <- D$U[[r]][, w.t, drop=FALSE]
              U.1 <- U.t[, w.1, drop = FALSE]
              Omega.1 <- diag(p.1) + K[r, r] * t(U.1) %*% U.1
              lambda.1 <- K[r, r] * t(Y.residual) %*% U.1 %*% solve(Omega.1)
              post.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1) ## Integrated probability
              ##-----------------------------------------------------------------------
            }
            ##------- Now update gamma -----------
            if(log(runif(1)) < post.1 + log(D$Theories.prob[[r]][t]) - log(1 - D$Theories.prob[[r]][t]) )
            {
              theta$gamma[[r]][t] <- 1
              w.r.t <- w.t[which(M.curr == 1)]
              theta$beta[[r]][w.r.t] <- rmvnorm.precision(lambda.1,Omega.1)
              w.r <- c(w.r,w.r.t)
            }else{
              theta$gamma[[r]][t] <- 0
              theta$beta[[r]][w.t] <- 0
            }
            ##------------------------------------
          }
        }
      }
    ##------- End looping through non-intercept theories ----
    
    ##---- Make a joint update of the regression parameter ---
    Y.residual <- sur.bma.theory.residual(theta,D,r)
    p.r <- length(w.r)
    U.r <- D$U[[r]][, w.r, drop = FALSE]
    Omega.r <- diag(p.r) + K[r, r] * t(U.r) %*% U.r
    lambda.r <- K[r, r] * t(Y.residual) %*% U.r %*% solve(Omega.r)
    theta$beta[[r]] <- rep(0, dim(D$U[[r]])[2])
    theta$beta[[r]][w.r] <- rmvnorm.precision(lambda.r, Omega.r)
    ##--------------------------------------------------------
    
    ##----- Retrieve residual --------------
    eps[,r] <- D$Y[,r] - D$U[[r]] %*% theta$beta[[r]]
    ##--------------------------------------
  }
  ##------- End Looping through equations ------------

  ##-------- Resample precision -------------
  E <- (D$n - 1) * cov(eps)
  theta$K <- rwish(D$n + 3, diag(D$R) + E)
  ##-----------------------------------------
  
  return(theta)
}

sur.bma.theory.init <- function(D)
  {

    ##------- Init Theta ------
    theta <- NULL
    theta$gamma <- list()
    theta$M <- list()
    theta$beta <- list()
    ##--------------------------
    
    ##-------- Get residuals --------
    eps <- matrix(0, D$n, D$R)
    ##-------------------------------

    ##--- Decide whether theories are included ------
    for(r in 1:D$R)
    {
      n.t <- D$n.t[r]
      theta$M[[r]] <- list()
      theta$gamma[[r]] <- rbinom(n.t,1,D$Theories.prob[[r]])
      theta$gamma[[r]][D$which.intercept[r] ] <- 1
      p.r <- dim(D$U[[r]])[2]
      theta$beta[[r]] <- rep(0,p.r)
      w.r <- NULL

      ##-- Intercept theory is the "always include" set so treat differently ---
      w.t <- D$w.t[[r]][[ D$which.intercept[r] ]]
      theta$M[[r]][[ D$which.intercept[r] ]] <- rep(1,length(w.t))
      w.r <- c(w.r, w.t)
      ##---------------------------------------------------

      ##---- Now set up the rest of the theories -----
      for(t in 1:n.t)
      {
        if(! (t == D$which.intercept[r]) )
        {
          w.t <- D$w.t[[r]][[t]]
          theta$M[[r]][[t]] <- rbinom(length(w.t),1,D$pi.M[[r]][w.t])
          if(all(theta$M[[r]][[t]] == 0))
          {
            theta$M[[r]][[t]] <- c(1,rep(0,length(w.t) - 1))
          }
          w.M <- which(theta$M[[r]][[t]] == 1)
          if(theta$gamma[[r]][t])
          {
            w.r <- c(w.r, w.t[w.M]) ## I just vommitted a little.
          }
        }
      }
      ##----------------------------------------------

      ##-------- Now sample model level beta --------
      U.r <- D$U[[r]][,w.r,drop=FALSE]
      p.r <- length(w.r)
      Y <- D$Y[,r]
      Omega.r <- diag(p.r) + t(U.r) %*% U.r
      lambda.r <- t(Y) %*% U.r %*% solve(Omega.r)
      theta$beta[[r]] <- rep(0, dim(D$U[[r]])[2])
      theta$beta[[r]][w.r] <- rmvnorm.precision(lambda.r, Omega.r)
      eps[,r] <- Y - D$U[[r]] %*% theta$beta[[r]]
      ##-----------------------------------------------
    }
    ##-------- End loop through equations ------

    ##-------- Resample precision -------------
    E <- (D$n - 1) * cov(eps)
    theta$K <- rwish(D$n + 3, diag(D$R) + E)
    ##-----------------------------------------
    
    return(theta)
  }

sur.bma.theory.results.init <- function(D,odens)
  {
    
    R <- D$R
    n <- D$n

    ##-------- Information to be returned ----------
    results <- NULL
    results$beta <- list()
p    results$beta.bar <- list()
    results$M <- list()
    results$M.bar <- list()
    results$gamma <- list()
    results$gamma.bar <- list()
    for(r in 1:R)
      {
        p.U.r <- dim(D$U[[r]])[2]
        results$beta[[r]] <- matrix(0, odens, p.U.r)
        results$beta.bar[[r]]  <- rep(0,p.U.r)
        results$M[[r]] <- matrix(0, odens, p.U.r)
        results$M.bar[[r]]  <- rep(0,p.U.r)
        results$gamma[[r]] <- matrix(0, odens,D$n.t[r])
        results$gamma.bar[[r]] <- rep(0, D$n.t[r])
      }
    results$K <- array(dim = c(R,R,odens))
    results$K.bar <- matrix(0,R,R)
    ##--------------------------------------------
    
    return(results)
  }

sur.bma.theory <- function(Y, U, s=1e3, b = round(s/10),
                           full = FALSE,odens = min(c(5e3,s-b)),
                           pi.M = NULL,
                           Theories = NULL,
                           Theories.prob = NULL,
                           print.every = round(s/10),
                           which.intercept = 1)
  {

    ##----- Setup data object ----------
    D <- NULL
    D$Y <- as.matrix(Y)
    D$n <- dim(D$Y)[1]
    D$R <- length(U)
    D$U <- list()
    R <- D$R
    D$p.r <- NULL
    for(r in 1:R)
    {
      D$U[[r]] <- as.matrix(U[[r]])
      D$p.r[r] <- dim(D$U[[r]])[2]
    }
    if(length(which.intercept) == 1)
    {
      D$which.intercept <- rep(which.intercept, D$R)
    }else{
      D$which.intercept <- which.intercept
    }
    ##---------------------------------


    ##--------- Setup Theories --------
    C <- list()
    D$n.t <- NULL
    D$w.t <- list()
    for(r in 1:R)
    {
      C[[r]] <- list()
      n.theory <- length(unique(Theories[[r]]))
      D$n.t[r] <- n.theory
      D$w.t[[r]] <- list()
      for(t in 1:n.theory)
      {
        w.theory <- which(Theories[[r]] == t)
        D$w.t[[r]][[t]] <- w.theory
        if(t != D$which.intercept[r] )
        {
          U.theory <- U[[r]][,w.theory,drop=FALSE]
          C[[r]][[t]] <- cor(U.theory)
        }
      }
    }
    D$C <- C
    D$Theories <- Theories
    ##----------------------------------

    ##------ Setup Theories Prob -------
    if(is.null(Theories.prob))
    {
      Theories.prob <- list()
      for(r in 1:R)
      {
        n.t <- D$n.t[r]
        Theories.prob[[r]] <- rep(.5,n.t)
        Theories.prob[[r]][D$which.intercept[r]] <- 1
      }
    }
    D$Theories.prob <- Theories.prob
    ##----------------------------------

    ##----- Variable Level Theories ----
    if(is.null(pi.M))
    {
      pi.M <- list()
      for(r in 1:R)
      {
        pi.M[[r]] <- rep(0, D$p.r[r])
        for(t in 1:D$n.t[r])
        {
          w.t <- D$w.t[[r]][[t]]
          if(t == D$which.intercept[r])
            {
              pi.M[[r]][w.t] <- 1
            }else{
              pi.M[[r]][w.t] <- .5
            }
        }
      }
    }
    D$pi.M <- pi.M
    ##----------------------------------
    
    ##----- Setup --------------------
    theta <- sur.bma.theory.init(D)
    results <- sur.bma.theory.results.init(D,odens)
    ##---------------------------------

    ##------ Bookkeeping --------------
    which.save <- round(seq(b + 1, s, length = odens))
    save.loc <- 1
    next.save <- which.save[save.loc]
    ##-------------------------------

    ##------ Run! ------------------
    for(i in 1:s)
      {

        if(i %% print.every == 0)print(paste("Iteration", i,Sys.time()))
        theta <- sur.bma.theory.update(theta,D)
##        print(c(round(theta$nu, 3), round(var(theta$alpha),3), round(mean(theta$alpha),3)))
        ##---------- Record ----------------------------
        if(i == next.save)
          {
            for(r in 1:R)
              {
                w.r <- which(theta$beta[[r]] != 0)
                results$M[[r]][save.loc,w.r] <- 1
                results$beta[[r]][save.loc,] <- theta$beta[[r]]
                results$gamma[[r]][save.loc,] <- theta$gamma[[r]]
                results$K[,,save.loc] <- theta$K
              }
            save.loc <- save.loc + 1
            next.save <- which.save[save.loc]
          }
        if(i > b)
          {
            for(r in 1:R)
              {
                w.r <- which(theta$beta[[r]] != 0)
                temp <- rep(0, dim(D$U[[r]])[2])
                temp[w.r] <- 1
                results$M.bar[[r]] <- results$M.bar[[r]] + temp / (s - b)
                results$beta.bar[[r]] <- results$beta.bar[[r]] + theta$beta[[r]]/ (s - b)
                results$gamma.bar[[r]] <- results$gamma.bar[[r]] + theta$gamma[[r]]/ (s - b)
              }
            results$K.bar <- results$K.bar + theta$K / (s - b)
          }
        ##------------------------------------------------------
      }
    ##--------------------------------------------

    results$D <- D
    return(results)
  }
