rwish <- function (delta, D) 
{
  ##Bartlet Decomposition
  T <- chol(solve(D))
  p <- dim(D)[1]
  Psi <- matrix(0, p, p)
  for (i in 1:p) Psi[i, i] <- sqrt(rchisq(1, delta + p - i))
  if (p > 1)
    {
      for (i in 1:(p - 1))
        {
          for (j in (i + 1):p)
            {
              Psi[i, j] <- rnorm(1)
            }
        }
    }
  Psi <- Psi %*% T
  return(t(Psi) %*% Psi)
}

rmvnorm.precision <- function(mu,K)
  {
    p <- dim(K)[2]
    Q <- chol(K)
    z <- rnorm(p)
    x <- backsolve(Q,z) + mu
    return(x)
  }

mydet <- function(A)
  {
    return(2 * sum(log(diag(chol(A)))))
  }



sur.bma.update.r <- function(theta,D,r)
  {

    R <- D$R
    ##------ Stoopid ----------------
    l <- NULL
    l$M <- theta$M[[r]]
    l$beta <- theta$beta[[r]]
    Y <- D$Y[,r]
    U <- D$U[[r]]
    K <- theta$K
    for(i in 1:R)
      {
        if(i != r)
          {
            Y <- Y + K[i,r]/K[r,r] * (D$Y[,i] - D$U[[i]] %*% theta$beta[[i]])
          }
      }

    Y <- sqrt(theta$alpha) * Y
    
    ##------ Keep book --------
    n <- length(Y)
    p.U <- dim(U)[2]
    M <- l$M
    ##------------------------

    ##----- Flip a variable -----
    pp <- (1 - M) * D$pi.M[[r]] + M * (1 - D$pi.M[[r]])
    w <- sample(1:p.U,1,prob= pp)
    M.new <- M
    M.new[w] <- 1 - M.new[w]
    if(sum(M.new) == 0) return(l)
    ##-----------------------------
    
    ##---- New Score ------------
    p.2 <- sum(M.new)
    U.2 <- U[,(1:p.U)[M.new == 1],drop=FALSE]
    for(i in 1:n)
      {
        U.2[i,] <-sqrt(theta$alpha[i]) * U.2[i,]
      }
    Omega.2 <- diag(p.2) + K[r,r] * t(U.2) %*% U.2
    lambda.2 <- K[r,r] * t(Y) %*% as.matrix(U.2) %*% solve(Omega.2)
    score.2 <- 0.5 * lambda.2 %*% Omega.2 %*% t(lambda.2) - 0.5 * mydet(Omega.2)
    ##---------------------------

    ##---- Old Score ------------
    p.1 <- sum(M)
    U.1 <- U[,(1:p.U)[M == 1],drop=FALSE]
    for(i in 1:n)
      {
        U.1[i,] <- sqrt(theta$alpha[i]) * U.1[i,]
      }
    Omega.1 <- diag(p.1) + K[r,r] * t(U.1) %*% U.1
    lambda.1 <- K[r,r] * t(Y) %*% U.1 %*% solve(Omega.1)
    score.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1)
    ##---------------------------

    ##---- Decide ---------------
    alpha <- score.2 - score.1
    if(log(runif(1)) < alpha)
      {
        M <- M.new
        Omega.1 <- Omega.2
        lambda.1 <- lambda.2
      }
    ##---------------------------

    beta.1 <- rep(0,p.U)
    beta.1[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)
    l$M <- M
    l$beta <- beta.1

    return(l)
  }

sur.bma.update.alpha <- function(theta,D)
  {
    n <- D$n
    alpha <- rep(0,n)
    K <- theta$K
    eps <- matrix(0, D$n, D$R)
    R <- D$R
    for(r in 1:R)
      {
        eps[,r] <- D$Y[,r] - D$U[[r]] %*% theta$beta[[r]]
      }

    alpha[1] <- 1
    for(i in 2:n)
      {
        s <- eps[i,] %*% K %*% eps[i,]
        alpha[i] <- rgamma(1,(theta$nu + R)/2, (theta$nu + s)/2)
      }
    return(alpha)
  }

sur.bma.update.nu <- function(theta)
  {
    n <- length(theta$alpha)
    nu <- theta$nu
    alpha <- theta$alpha
    nu.prime <- rtnorm(1,nu,1,lower = 0)

    score.lik.top <- sum(dgamma(alpha,nu.prime/2,nu.prime/2,log=TRUE))
    score.prior.top <- dtnorm(nu.prime, 5,sqrt(100),lower = 0, log = TRUE)
    score.prop.top <- dtnorm(nu, nu.prime,1,lower = 0, log = TRUE)

    score.lik.bot <- sum(dgamma(alpha,nu/2,nu/2,log=TRUE))
    score.prior.bot <- dtnorm(nu, 5, sqrt(100), lower = 0, log = TRUE)
    score.prop.bot <- dtnorm(nu.prime, nu,1,lower = 0, log = TRUE)

    mh <- score.lik.top + score.prior.top + score.prop.top - score.lik.bot - score.prior.bot - score.prop.bot
    
    if(log(runif(1)) < mh)
      {
        return(nu.prime)
      }else{
        return(nu)
      }
  }

sur.bma.update <- function(theta, D,full)
  {
    n <- D$n
    eps <- matrix(0, D$n, D$R)
    R <- D$R
    for(r in 1:R)
      {
        l <- sur.bma.update.r(theta,D,r)
        theta$M[[r]] <- l$M
        theta$beta[[r]] <- l$beta
        eps[,r] <- D$Y[,r] - D$U[[r]] %*% theta$beta[[r]]
      }

    if(theta$dispersion)
      {
        for(i in 1:n)
          {
            eps[i,] <- sqrt(theta$alpha[i]) * eps[i,]
          }
      }

    U <- (D$n - 1) * cov(eps)
    theta$K <- rwish(D$n + 3, diag(R) + U)
    
    if(theta$dispersion)
      {
        theta$alpha <- sur.bma.update.alpha(theta,D)
        theta$nu <- sur.bma.update.nu(theta)
      }


    return(theta)
  }

sur.bma.init <- function(D,full,dispersion)
  {

    theta <- NULL
    theta$beta <- list()
    theta$M <- list()
    for(r in 1:D$R)
      {
        p.U.r <- dim(D$U[[r]])[2]
        theta$beta[[r]] <- solve(t(D$U[[r]]) %*% D$U[[r]]) %*% t(D$U[[r]]) %*% D$Y[,r]
        
        if(full)
          {
            theta$M[[r]] <- 1
          }else{
            theta$M[[r]] <- rbinom(p.U.r,1,D$pi.M[[r]])
          }
        theta$beta[[r]][theta$M[[r]] == 0] <- 0
      }
    theta$K <- diag(D$R)
    theta$G <- matrix(1,D$R,D$R)
    theta$alpha <- rep(1,D$n)
    theta$nu <- rnorm(1,20,.1)
    theta$dispersion <- dispersion
    return(theta)
  }

sur.bma.results.init <- function(D,odens)
  {

    R <- D$R
    n <- D$n
    ##-------- Information to be returned ----------
    results <- NULL
    results$beta <- list()
    results$beta.bar <- list()
    results$M <- list()
    results$M.bar <- list()
    for(r in 1:R)
      {
        p.U.r <- dim(D$U[[r]])[2]
        results$beta[[r]] <- matrix(0, odens, p.U.r)
        results$beta.bar[[r]]  <- rep(0,p.U.r)
        results$M[[r]] <- matrix(0, odens, p.U.r)
        results$M.bar[[r]]  <- rep(0,p.U.r)
      }
    results$K <- array(dim = c(R,R,odens))
    results$G <- array(dim = c(R,R,odens))
    results$K.bar <- matrix(0,R,R)
    results$G.bar <- matrix(0,R,R)
    results$nu <- rep(0, odens)
    results$nu.bar <- 0
    results$y.hat <- array(dim = c(n,odens,R))
    results$y.hat.bar <- matrix(0,n,R)
    ##--------------------------------------------
    
    return(results)
  }

sur.bma <- function(Y, U, s=1e3, b = round(s/10),
                    full = FALSE,odens = min(c(5e3,s-b)),
                    pi.M = NULL,
                    print.every = round(s/10),
                    dispersion = FALSE)
  {
    D <- NULL
    D$Y <- as.matrix(Y)
    D$n <- dim(D$Y)[1]
    D$R <- length(U)
    D$U <- list()
    R <- D$R
    for(r in 1:R) D$U[[r]] <- as.matrix(U[[r]])
    if(is.null(pi.M))
      {
        D$pi.M <- list()
        for(r in 1:R) D$pi.M[[r]] <- rep(.5,dim(U[[r]])[2])
      }else{
        D$pi.M <- pi.M
      }
    
    theta <- sur.bma.init(D,full,dispersion)
    results <- sur.bma.results.init(D,odens)
    which.save <- round(seq(b + 1, s, length = odens))
    save.loc <- 1
    next.save <- which.save[save.loc]

    for(i in 1:s)
      {

        if(i %% print.every == 0)print(i)
        theta <- sur.bma.update(theta,D)
##        print(c(round(theta$nu, 3), round(var(theta$alpha),3), round(mean(theta$alpha),3)))
        ##---------- Record ----------------------------
        if(i == next.save)
          {
            for(r in 1:R)
              {
                results$M[[r]][save.loc,] <- theta$M[[r]]
                results$beta[[r]][save.loc,] <- theta$beta[[r]]
                results$K[,,save.loc] <- theta$K
                results$G[,,save.loc] <- theta$G
                results$nu[save.loc] <- theta$nu
                results$y.hat[,save.loc,r] <- D$U[[r]] %*% theta$beta[[r]]
              }
            save.loc <- save.loc + 1
            next.save <- which.save[save.loc]
          }
        if(i > b)
          {
            for(r in 1:R)
              {
                results$M.bar[[r]] <- results$M.bar[[r]] + theta$M[[r]]/ (s - b)
                results$beta.bar[[r]] <- results$beta.bar[[r]] + theta$beta[[r]]/ (s - b)
                results$y.hat.bar[,r] <- results$y.hat.bar[,r] + (D$U[[r]] %*% theta$beta[[r]])/(s-b)
              }
            results$K.bar <- results$K.bar + theta$K / (s - b)
            results$G.bar <- results$G.bar + theta$G / (s - b)
            results$nu.bar <- results$nu.bar + theta$nu / (s- b)
          }
        ##------------------------------------------------------
      }
    results$D <- D
    return(results)
  }


sur.bma.summary <- function(a,U.lab)
  {

    tbl <- list()
    for(r in 1:length(a$beta))
      {
        tbl[[r]] <- cbind(a$M.bar[[r]], t(apply(a$beta[[r]],2,"quantile",c(.025,.5, .975))),a$beta.bar[[r]], apply(a$beta[[r]],2,sd))
        A <- NULL
        B <- NULL
        for(k in 1:dim(tbl[[r]])[1])
          {
            A[k] <- sd(a$beta[[r]][a$beta[[r]][,k] != 0,k])
            B[k] <- mean(a$beta[[r]][a$beta[[r]][,k] != 0,k])
          }
        tbl[[r]] <- cbind(tbl[[r]], cbind(B,A))
      }
    

                                        #changed rownames(tbl[[1]]) <- c(X.lab,"X1.sq","Intercept",W.lab)
    for(r in 1:length(U.lab))
      {
        rownames(tbl[[r]]) <- U.lab[[r]]
        colnames(tbl[[r]]) <- c("Prob","Low","Med","Up","Mean","SD","CondMean","CondSD")
      }

    return(tbl)
  }

## a is an object returned from sur.bma
## I assume the first element in Y and U is the dependent equation
## which.z is a list of which objects in the remaining covariates are instruments
sur.bma.iv.diagnostics <- function(a, Z)
  {
    Z <- as.matrix(Z)
    results <- NULL
    results$BTIV.collapsed <- NULL
    results$BTIV.collapsed.cond <- NULL
    q <- dim(Z)[2]
    n.0 <- 1
    n <- dim(a$D$Y)[1]
    R <- dim(a$D$Y)[2]
    eps <- matrix(0,n,R)
    
    for(i in 1:dim(a$beta[[1]])[1])
        {
          K <- a$K[,,i]
          for(r in 1:dim(a$D$Y)[2])eps[,r] <- (a$D$Y[,r] - a$D$U[[r]] %*% a$beta[[r]][i,])
          eps.hat <- eps[,1]
          eps.cond <- eps[,1] + eps[,2:R] %*% K[1,2:R]/K[1,1]
          K.Z <- (n.0 * diag(q) + t(Z) %*% Z)
          beta.Z <- solve(K.Z) %*% t(Z) %*% eps.hat
          S.Z <- t(eps.hat) %*% eps.hat - t(beta.Z) %*% K.Z %*% beta.Z
          score.top <- (q/2) * log(n.0) - 1/2 * mydet(K.Z) - (1 + n)/2 * log((1/2 + S.Z/2))
          score.bottom <- -(1 + n)/2 * log((1/2 + t(eps.hat) %*% eps.hat/2))
          results$BTIV.collapsed[i] <- 1/(1 + exp(score.bottom - score.top))

          beta.Z <- solve(K.Z) %*% t(Z) %*% eps.cond
          S.Z <- t(eps.cond) %*% eps.cond - t(beta.Z) %*% K.Z %*% beta.Z
          score.top <- (q/2) * log(n.0) - 1/2 * mydet(K.Z) - (1 + n)/2 * log((1/2 + S.Z/2))
          score.bottom <- -(1 + n)/2 * log((1/2 + t(eps.cond) %*% eps.cond/2))
          results$BTIV.collapsed.cond[i] <- 1/(1 + exp(score.bottom - score.top))
        }
        return(results)

  }

sur.bma.poisson.init <- function(D,full,dispersion)
    {
        theta <- NULL
        theta$D.normal <- NULL
        theta$D.normal$U <- D$U
        theta$D.normal$Y <- D$Y
        theta$D.normal$Y[,1] <- log(D$Y[,1])
        theta$D.normal$Y[which(D$Y[,1] == 0),1] <- log(.5)
        theta$D.normal$R <- D$R
        theta$D.normal$n <- D$n
        theta$D.normal$pi.M <- D$pi.M
        theta$theta.normal <- sur.bma.init(theta$D.normal,full, dispersion)
        return(theta)

    }

sur.bma.poisson.results.init <- function(D, odens)
    {
        results <- NULL
        results$results.normal <- sur.bma.results.init(D,odens)
        results$RE <- matrix(0, odens,dim(D$Y)[1])
        results$updated <- rep(0, D$n)
        return(results)
    }

sur.bma.poisson.update.RE <- function(theta,D)
    {
        updated <- rep(0, D$n)
        beta.Y <- theta$theta.normal$beta[[1]]
        U <- theta$D.normal$U[[1]]
        mu.hat <- U %*% beta.Y

        epsilon <- matrix(0,D$n,D$R)
        for(r in 1:R)
            {
                beta <- theta$theta.normal$beta[[r]]
                U <- theta$D.normal$U[[r]]
                mu <- U %*% beta
                epsilon[,r] <- theta$D.normal$Y[,r] - mu
            }
        
        Y <- D$Y[,1]
        for(i in 1:n)
            {
                e.i <- epsilon[i,1]
                mu.i <- mu.hat[i]
                kappa <- theta$theta.normal$K[1,1]
                eta.i <- sum(-theta$theta.norma$K[1,-1]/kappa * epsilon[i,-1])
                f.prime <- -exp(mu.hat[i] + e.i) + Y[i] - kappa * (e.i - eta.i)
                f.2.prime <- -exp(mu.hat[i] + e.i) -kappa
                b <- f.prime - f.2.prime * e.i
                c <- -f.2.prime
                e.new <- rnorm(1, b/c, sqrt(1/c))

                f.prime.new <- -exp(mu.hat[i] + e.new) + Y[i] - kappa * (e.new - eta.i)
                f.2.prime.new <- -exp(mu.hat[i] + e.new) -kappa
                b.new <- f.prime.new - f.2.prime.new * e.new
                c.new <- -f.2.prime.new

                l.top <- dpois(Y[i],exp(mu.i + e.new), log=TRUE) + dnorm(e.new, eta.i, sqrt(1/kappa), log=TRUE) + dnorm(e.i, b.new/c.new, sqrt(1/c.new), log=TRUE)
                l.bottom <- dpois(Y[i],exp(mu.i + e.i), log=TRUE) + dnorm(e.i, eta.i, sqrt(1/kappa), log=TRUE) + dnorm(e.new, b/c, sqrt(1/c), log=TRUE)
                if(log(runif(1)) < l.top - l.bottom)
                    {
                        epsilon[i,1] <- e.new
                        updated[i] <- 1
                    }
            }
        theta$D.normal$Y[,1] <- mu.hat + epsilon[,1]
        theta$updated <- updated
        return(theta)
    }

sur.bma.poisson.update <- function(theta, D)
    {
        theta$theta.normal <- sur.bma.update(theta$theta.normal, theta$D.normal,theta$full)
        theta <- sur.bma.poisson.update.RE(theta, D)
        return(theta)
    }

sur.bma.poisson <- function(Y,U,s=1e3,b=round(s/10),full = FALSE,odens = min(c(5e3,s-b)),
                    pi.M = NULL,
                    print.every = round(s/10),
                    dispersion = FALSE)
    {

        D <- NULL
        D$Y <- as.matrix(Y)
        D$n <- dim(D$Y)[1]
        D$R <- length(U)
        D$U <- list()
        R <- D$R
        for(r in 1:R) D$U[[r]] <- as.matrix(U[[r]])
        if(is.null(pi.M))
            {
                D$pi.M <- list()
                for(r in 1:R) D$pi.M[[r]] <- rep(.5,dim(U[[r]])[2])
            }else{
                D$pi.M <- pi.M
            }
        
        theta <- sur.bma.poisson.init(D,full,dispersion)
        results <- sur.bma.poisson.results.init(D,odens)
        
        which.save <- round(seq(b + 1, s, length = odens))
        save.loc <- 1
        next.save <- which.save[save.loc]
        
        for(i in 1:s)
            {
                
                if(i %% print.every == 0)print(i)
                theta <- sur.bma.poisson.update(theta,D)
                ##        print(c(round(theta$nu, 3), round(var(theta$alpha),3), round(mean(theta$alpha),3)))
                ##---------- Record ----------------------------
                if(i == next.save)
                    {
                        for(r in 1:R)
                            {
                                results$results.normal$M[[r]][save.loc,] <- theta$theta.normal$M[[r]]
                                results$results.normal$beta[[r]][save.loc,] <- theta$theta.normal$beta[[r]]
                                results$results.normal$K[,,save.loc] <- theta$theta.normal$K
                                results$results.normal$G[,,save.loc] <- theta$theta.normal$G
                                results$results.normal$nu[save.loc] <- theta$theta.normal$nu
                                results$results.normal$y.hat[,save.loc,r] <- D$U[[r]] %*% theta$theta.normal$beta[[r]]
                            }
                        results$RE[save.loc,] <- theta$D.normal$Y[,1] - theta$D.normal$U[[1]] %*% theta$theta.normal$beta[[1]]
                        save.loc <- save.loc + 1
                        next.save <- which.save[save.loc]
                    }
                if(i > b)
                    {
                        for(r in 1:R)
                            {
                                results$results.normal$M.bar[[r]] <- results$results.normal$M.bar[[r]] + theta$theta.normal$M[[r]]/ (s - b)
                                results$results.normal$beta.bar[[r]] <- results$results.normal$beta.bar[[r]] + theta$theta.normal$beta[[r]]/ (s - b)
                                results$results.normal$y.hat.bar[,r] <- results$results.normal$y.hat.bar[,r] + (theta$D.normal$U[[r]] %*% theta$theta.normal$beta[[r]])/(s-b)
                            }
                        results$K.bar <- results$K.bar + theta$K / (s - b)
                        results$G.bar <- results$G.bar + theta$G / (s - b)
                        results$nu.bar <- results$nu.bar + theta$nu / (s- b)
                        results$updated <- results$updated + theta$updated / (s-b)
                        
                    }
                ##------------------------------------------------------
            }
        results$D <- D
        return(results)
    }
