summary.ivbma <- function(object,nms.U=NULL,nms.V=NULL,...)
  {
    x <- object
    p.U <- dim(x$lambda.bar)[1]
    r <- dim(x$lambda.bar)[2]
    p.V <- length(x$rho.bar)
    tbl.s1 <- list()
    tbl.s2 <- matrix(0,p.V,5)

    
    for(j in 1:r)
      {
        tbl.s1[[j]] <- matrix(0,p.U,5)
        for(i in 1:p.U)
          {
            tbl.s1[[j]][i,] <- c(x$M.bar[i,j],x$lambda.bar[i,j],quantile(x$lambda[i,j,],c(.025,.5,.975)))
          }
        if ( !is.null (nms.U))
          {
            rownames(tbl.s1[[j]]) <- nms.U
          }
        colnames(tbl.s1[[j]]) <- c("Prob","Mean", "Lower","Med","Upper")
      }


    for(i in 1:p.V)
      {
        tbl.s2[i,] <- c(x$L.bar[i], x$rho.bar[i],quantile(x$rho[,i],c(.025,.5,.975)))
      }
    if ( !is.null (nms.V))
      {
        rownames(tbl.s2) <- nms.V
      }
    colnames(tbl.s2) <- c("Prob","Mean","Lower","Med","Upper")

    l <- NULL
    l$tbl.s1 <- tbl.s1
    l$tbl.s2 <- tbl.s2
    l$tbl.cov <- x$Sigma.bar
    if(x$run.diagnostics)
      {
        l$Bayes.Sargan <- x$Bayesian.Sargan
        l$Sargan <- x$Sargan
      }
    return(l)
  }



collapse.model <- function(m)
  {
    return(paste(m,collapse=""))
  }

parse.model <- function(mod.str,nms)
  {
    p <- nchar(mod.str)
    mod <- ""
    k <- 1
    for(i in 1:p) ## Vectorize
      {
        tmp <- as.numeric(substr(mod.str,i,i))
        if(tmp == 1)
          {
            if(k == 1)
              {
                mod <- nms[i]
              }else{
                mod <- paste(mod,nms[i],sep=",")
              }
            k <- k + 1
          }
      }
    return(mod)
  }

model.tabulation <- function(M,nms,n=5)
  { 

    M.collapse <- apply(M,1,"collapse.model")
    table.collapse <- table(M.collapse) / sum(table(M.collapse))
    o.collapse <- order(table.collapse, decreasing=TRUE)
    if(length(o.collapse) > n)
    {
      top.n <- table.collapse[o.collapse[1:n]]
    }else{
      top.n <- table.collapse[o.collapse]
    }
    nms.top <- names(top.n)
    mod.top <- unlist(lapply(nms.top,"parse.model",nms))
    p.top <- unname(top.n)
    return(list(models = mod.top, prob = p.top))
  }
