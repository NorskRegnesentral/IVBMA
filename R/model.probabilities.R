model.priors.from.theories <- function(Theories, Theories.prob, U)
  {

    ##-------- Error checking -----------
    R <- length(Theories)
    if( (length(Theories.prob) != R)  | (length(U) != R) )
      {
        stop(print("The length of the Theories, Theories.prob and U objects are not all identical"))
      }
    p <- rep(NA,R)
    for(i in 1:R)
      {
        p[i] <- length(Theories[[i]])
        if(p[i] != dim(U[[i]])[2])
          {
            stop(print(paste("Theories i has",p[i],"elements but U[[i]] has",dim(U[[i]])[2],
                             "variables.  These need to be equal")))
          }
      }
    d <- rep(NA,R)
    for(i in 1:R)
      {
        d[i] <- length(Theories.prob[[i]])
        if(d[i] != max(Theories[[i]]))
          {
            stop(print(paste("Theories i has",max(Theories[[i]]),"theories in it, but Theories.prob[[i]] has",d[i],"elements.  These need to be equal")))
          }
      }
    ##------------------------------------

    ##----- Now backout the pi.M object ----
    pi.M <- list()
    for(i in 1:R)
      {
        pi.M[[i]] <- rep(NA,p[i])
        for(j in 1:d[i])
          {
            w.j <- which(Theories[[i]] == j)
            if(length(w.j) > 0)
              {
                pi.M[[i]][w.j] <- 1 - (1 - Theories.prob[[i]][j])^(1 / length(w.j))
              }
          }
      }
    ##---------------------------------------

    return(pi.M)
  }

posterior.theory.probability <- function(M,Theories)
  {
    R <- length(Theories)
    prob.theory <- list()
    for(i in 1:R)
      {
        d <- max(Theories[[i]])
        prob.theory[[i]] <- rep(NA,d)
        for(j in 1:d)
          {
            w.j <- which(Theories[[i]] == j)
            if(length(w.j) > 0)
              {
                M.j <- M[[i]][,w.j, drop=FALSE]
                M.in <- (rowSums(M.j) > 0) * 1
                prob.theory[[i]][j] <- mean(M.in)
              }
          }
      }

    return(prob.theory)
  }
