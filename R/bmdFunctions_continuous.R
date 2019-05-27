#' Main function for calculations with continuous response type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all 
#' @export
f.con <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.con")
  
  ans.all <- f.start.con(ans.all = ans.all, fitted = FALSE, track = track)
  
  # Change start values of model parameters
  if(!is.null(ans.all$startValues))
    ans.all$main.ans <- c(3, ans.all$main.ans)
  
  
  for(main.ans.single in c(ans.all$main.ans, 13)){
    
    switch(as.character(main.ans.single), 
        
        "3" = { # Change start values of model parameters
          
          ans.dd <- 2
#          if (max(ans.all$fct5) > 1) ans.dd <- menu(c("yes", 
#                    "no"), title = "\nDo you want to use previous MLE from model with one d as start value?\n")
#          
#          if (ans.dd == 1) {
#            
#            nr.dd <- max(ans.all$fct5)
#            dd.start <- ans.all$MLE[length(ans.all$MLE)]
#            ans.all$par.start <- c(ans.all$MLE, rep(dd.start, 
#                    nr.dd - 1))
#            ans.all$scale.dum <- c(ans.all$scale.dum, rep(1, 
#                    nr.dd - 1))
#            
#          }
          
          ans.all <- f.start.con(ans.all, adjust = TRUE, fitted = FALSE, track = track)
          list.logic <- F
          
        },
        
        "4" = { # Fit model
          
          ans.all <- f.mm4.con(ans.all, track = track)
          
        },
        
        "6" = { # Calculate CED
          
          if (!ans.all$fitted) {
            
            stop("First fit the model")
            
          }
          
          if (ans.all$model.ans %in% c(1, 11)) {
            
            stop("BMD not defined for null or full model")
            
          }
          
          ans.all <- f.mm6.con(ans.all, track = track)
          
        },
        
        "7" = { # Calculate CED conf interval with bootstrap
          
          if (is.na(ans.all$CED[1])) cat("\nYou did not calculate CED! First use option 6 from main menu\n") else {
            ans.all$plot.ans <- 1
            ans.all <- f.mm7.con(ans.all, track = track)
          }
          
        },
        
        "13" = { # Return results
          
#					main.ans <- 13
          if (track) 
            print("f.con:  END")
          return(ans.all)
          
        }) 
    
  }
  
}



#' Define initial parameter values for the chosen model, continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param adjust boolean TRUE if start values of the parameters should be adjusted;
#' default value is FALSE
#' @param fitted boolean TRUE if the model has already been fitted;
#'  default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.start.con <- function (ans.all, adjust = FALSE, fitted = FALSE, track = FALSE) {
  
  if (track) 
    print("f.start.con")
  
  cc.upp <- 10000
  dd.upp <- 10
  
  ans.all$fitted <- fitted
  ans.all$nr.aa <- max(ans.all$fct1)
  ans.all$nr.bb <- max(ans.all$fct2)
  ans.all$nr.var <- max(as.numeric(ans.all$fct3))
  ans.all$nr.cc <- max(ans.all$fct4)
  ans.all$nr.dd <- max(ans.all$fct5)
  
  with(ans.all, {
        
        if (dtype == 5) dtype <- 1
        
        
        y.sort <- yy[order(x)]
        x.sort <- sort(x)
        
        if (!adjust) {
          
          if (dtype %in% c(10, 15, 250, 260) & length(x) < 10) {
            # Continuous summary data
            
            xtmp <- x.sort
            ytmp <- y.sort
            
          } else {
            # Make 5 groups for linear model fit
            
            sub <- round(length(x.sort)/5)
            if (sub == 0) 
              sub <- 1
            xtmp <- numeric()
            ytmp <- numeric()
            e1 <- 0
            e2 <- 0
            for (k in 1:5) {
              e2 <- e2 + sub
              qq <- x.sort[e1:e2]
              xtmp[k] <- mean(qq, na.rm = TRUE)
              qq <- y.sort[e1:e2]
              if (dtype %in% c(25, 250)) 
                ytmp[k] <- mean(qq, na.rm = TRUE)
              if (dtype %in% c(26, 260)) 
                ytmp[k] <- mean(sqrt(qq, na.rm = TRUE))^2
              else if (!(dtype %in% c(25, 26, 250, 260))) 
                ytmp[k] <- exp(mean(log(qq[!is.na(qq)])))
              e1 <- e1 + sub
            }
          }
          if (length(ytmp[!is.na(ytmp)]) == 3) {
            ytmp[4:5] <- ytmp[3]
            xtmp[4:5] <- xtmp[3]
          }
          if (length(ytmp[!is.na(ytmp)]) == 4) {
            ytmp[5] <- ytmp[4]
            xtmp[5] <- xtmp[4]
          }
          top <- length(ytmp)
          aa <- NA
          bb <- NA
          cc <- NA
          dd <- NA
          
          if (model.ans != 1) {
            
            lm.res <- f.start.lm.con(xtmp, ytmp, model.ans, dtype, track = track)
            
            aa <- abs(lm.res$intercept)
            if (model.ans %in% c(2, 3, 7, 8, 17, 18)){
              
              bb <- lm.res$slope
              
            } else {
              
              bb <- abs(lm.res$slope)
              
            }
            
            increase <- lm.res$increase
            
          } else {
            
            increase <- 1
            
          }
          
          #ans.all$CES <- CES*increase
          CES <- CES*increase
          
          switch(as.character(model.ans), 
              
              "1" = { #Null model
                
                aa <- mean(ytmp)
                regr.par <- rep(aa, nr.aa)
                
              }, 
              
              "11" = { #Full model
                
                if (dtype %in% c(10, 15, 250, 260)) {
                  
                  if (dtype %in% c(10, 15)) regr.par <- exp(mn.log)
                  if (dtype == 250) regr.par <- (mn.log)
                  if (dtype == 260) regr.par <- (mn.log)^2
                  
                } else {
                  
                  regr.par <- numeric()
                  nn <- numeric()
                  if (nr.var == 1) { 
                    
                    for (jj in levels(factor(fct2))) for (ii in levels(factor(fct1))) {
                        x.tmp <- x[fct1 == ii & fct2 == jj]
                        y.tmp <- yy[fct1 == ii & fct2 == jj]
                        regr.tmp <- as.numeric(exp(tapply(log(y.tmp), 
                                    x.tmp, mean)))
                        regr.par <- c(regr.par, regr.tmp)
                        nn.tmp <- tapply(y.tmp, x.tmp, length)
                        nn <- c(nn, nn.tmp)
                      }
                    
                  } else {
                    
                    regr.par <- numeric()
                    for (jj in levels(factor(fct3))) {
                      x.tmp <- x[fct3 == jj]
                      y.tmp <- yy[fct3 == jj]
                      regr.tmp <- as.numeric(exp(tapply(log(y.tmp), 
                                  x.tmp, mean)))
                      regr.par <- c(regr.par, regr.tmp)
                      nn.tmp <- tapply(y.tmp, x.tmp, length)
                      nn <- c(nn, nn.tmp)
                    }
                  }
                  ans.all$nn <- nn
                }
                
                ans.all$lower <- regr.par
                ans.all$upper <- regr.par
                
              }, 
              
              "13" = { #E3 - CED
                
                if (fitted) {
                  
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  
                }
                
                if (!fitted || is.na(bb) || bb == 0) {
                  
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  CED <- (1/bb) * log(CES + 1)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                          min(x))/3
                }
                dd <- 1
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                    rep(dd, nr.dd))
                
              }, 
              
              "14" = { #E4 - CED (used as full model for dtype = 3)
                
                if (fitted) {
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  cc <- 0
                }
                if (!fitted || is.na(bb) || bb == 0) {
                  cc <- ytmp[top]/aa
                  if (CES < 0 & cc > 1 + CES) {
                    cc <- 1 + CES - (1 + CES)/100
                  }
                  if (CES > 0 & cc < 1 + CES) {
                    cc <- 1 + CES + CES/100
                  }
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  CED <- -(1/bb) * logb((CES + 1 - cc)/(1 - cc))
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                          min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                    rep(cc, nr.cc))
                if (cc < 1) {
                }
                if (cc > 1) {
                }
                
              },
              
              "15" = { #E5 - CED
                
                if (fitted) {
                  
                  aa <- par.start[nr.var + 1]
                  bb <- par.start[nr.var + nr.aa + 1]
                  
                }
                if (!fitted || is.na(bb) || bb == 0) {
                  
                  dd <- 1
                  cc <- ytmp[top]/aa
                  if ((CES < 0) & (cc > 1 + CES)) {
                    cc <- 1 + CES - (1 + CES)/100
                  }
                  if ((CES > 0) & (cc < 1 + CES)) {
                    cc <- 1 + CES + CES/100
                  }
                  if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                  if (cc == 0) cc <- 0.1
                  if (cc == 1) cc <- 0.9
                  CED <- (-(1/bb) * log((-CES + 1 - cc)/(1 - 
                                  cc)))^(1/dd)
                  CED <- abs(CED)
                  if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                          min(x))/3
                }
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                    rep(cc, nr.cc), rep(dd, nr.dd))
              }, 
              
              "23" = { #H3 - CED
                if (CES > 0) bb <- -bb
                dd <- 1
                CED <- f.inv.con(model.ans = 18, c(aa, bb, dd), CES, track = track)
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                    rep(dd, nr.dd))
                y.expect <- f.expect.con(model.ans, xtmp, regr.par, 
                    CES = CES, increase = increase, track = track)
                test <- abs(sum(sign(y.expect)))
                if (!any(is.na(test))) if (test != length(y.expect)) {
                    bb <- 1.5 * bb
                    regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(dd, nr.dd))
                  }
                
              }, 
              
              "25" = { # H5 - CED
                
                cc <- ytmp[top]/aa
                dd <- 1
                if (CES < 0 & cc > 1 + CES) {
                  cc <- 1 + CES - (1 + CES)/100
                }
                if (CES > 0 & cc < 1 + CES) {
                  cc <- 1 + CES + CES/100
                }
                if (is.na(bb) || bb == 0) bb <- mean(xtmp)
                if (cc == 0) cc <- 0.1
                if (cc == 1) cc <- 0.9
                CED <- f.inv.con(model.ans = 20, c(aa, bb, cc, dd), CES,
                    track = track)
                if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                        min(x))/3
                regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                    rep(cc, nr.cc), rep(dd, nr.dd))
              })
          
          ans.all$regr.par <- regr.par
          text.par <- f.text.par(ans.all)
          
          if (!cont) {
            if (track) print("f.start.con:  END.sub")
            
            return(ans.all)
          }
          
          ans.all$nrp <- length(text.par) - nr.var
          
          if (!(dtype %in% c(10, 15, 250, 260))) {
            
            expect.trans <- f.expect.con(model.ans, x, regr.par, 
                fct1 = fct1, fct2 = fct2, fct3 = fct3, CES = CES, 
                twice = twice, ttt = 0, y = yy, increase = increase, 
                x.mn = x.mn, par.start = par.start, track = track)
            
            yy.trans <- yy
            
            if (dtype == 26) {
              
              expect.trans <- sqrt(expect.trans)
              yy.trans <- sqrt(yy)
              
            } else if (dtype %in% c(1, 5, 15)) {
              
              expect.trans <- log(expect.trans)
              yy.trans <- log(yy)
              
            }
            
            if (!model.ans %in% c(1:10, 12:20, 22:25)) {
              
              ans.all$lower <- c(1e-06, ans.all$lower)
              ans.all$upper <- c(1, ans.all$upper)
              
            }
            
            if (any(is.na(expect.trans))) {
              
              ans.all$par.start <- c(rep(NA, nr.var), regr.par)
              warning("Adjust start values")
              ans.all$adjust.start <- TRUE
              
            }
            
            var.start <- var(yy.trans - expect.trans)
            
          } else {
            
            var.start <- mean(sd2.log)
            
          }
          
          ans.all$par.start <- c(rep(var.start, nr.var), regr.par)
          
          loglik.first <- -f.lik.con(ans.all$par.start, 
              x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
              sd2.log, nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
              ttt = 0, fct4 = fct4, fct5 = fct5, 
              cens.up = cens.up, par.tmp = NA, increase = increase, 
              x.mn = x.mn, track = track)
          
          if (is.na(loglik.first) | loglik.first < -1e+05) {
            
            warning("Adjust start values before fitting the model")
            ans.all$adjust.start <- TRUE
            return(ans.all)
            
          }
          
          ans.all$loglik.first <- loglik.first
          
        } else { # if adjust == TRUE
          
          loglik.new <- NA
          ans.all$par.start <- startValues
          
          if (sum(is.finite(par.start)) == length(par.start)) 
            loglik.new <- -f.lik.con(par.start, 
                x, y, dtype, fct1, fct2, fct3, model.ans, 
                mn.log, sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                twice = twice, ttt = ttt, fct4 = fct4, fct5 = fct5, 
                cens.up = cens.up, par.tmp = NA, increase = increase, 
                x.mn = x.mn, track = track)
          
          if (model.ans == 11 & track) 
            cat("\nATTENTION: log-likelihood not implemented for full model\n")
          
          if(is.na(loglik.new))
            ans.all$errorAdjustStartValues <- TRUE
#            stop("Please choose other start values for the model parameters:\n The log-likelihood could not be calculated")
          
          ans.all$regr.par <- par.start[-(1:nr.var)]
          ans.all$loglik.old <- loglik.new          
          ans.all$adjust.start <- FALSE
          ans.all$loglik <- loglik.new
          
        }
        
        ans.all$nr.var <- nr.var
        ans.all$npar <- length(ans.all$par.start)
        ans.all$CED <- NA
        ans.all$increase <- increase
        ans.all$CES <- CES
        ans.all$max.lev <- nr.aa * nr.bb
        
        ans.all <- f.constr.con(ans.all, track = track)
        
        if (track)	print("f.start.con:  END")
        
        return(ans.all)
        
      })
}



#' Calculate likelihood for continuous response
#' @param theta numeric vector, the initial regression parameter values
#' @param x numeric vector, the dose values
#' @param y numeric vector, the response values 
#' @param dtype integer, determines the type of response
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param fct3 numeric, value for parameter var
#' @param model.ans integer, determines the model that will be fitted 
#' @param mn.log numeric vector, transformation of the response values,
#' see f.execute()
#' @param sd2.log  numeric vector, transformation of the sd of the response
#' values, see f.execute()
#' @param nn numeric vector, the number of responses per dose level, for 
#' continuous summary data
#' @param Vdetlim numeric vector, values of detection limit
#' @param CES numeric, value for the CES
#' @param twice boolean, if some parameter values are equal, see f.execute()
#' @param ttt numeric, time variable 
#' @param fct4 numeric, value for parameter c
#' @param fct5 numeric, value for parameter d
#' @param cens.up numeric, value for right censoring
#' @param lb numeric vector, determines the lower bound for theta;
#' default value is -Inf
#' @param ub numeric vector, determines the upper bound for theta;
#' default value is Inf
#' @param par.tmp numeric vector, regression parameter values, see f.pars() 
#' @param increase boolean, whether the response values are increasing or 
#' decreasing for increasing dose values 
#' @param x.mn numeric value, the mean of dose values, see f.execute()
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric value, minus the sum of the scores (total likelihood)
#' @export
f.lik.con <- function (theta, x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
    sd2.log, nn, Vdetlim, CES, twice = TRUE, ttt = 0,  
    fct4 = 1, fct5 = 1, cens.up = NA, lb = -Inf, ub = Inf, par.tmp, 
    increase = increase, x.mn = NA, track = FALSE) {
  
  if (track) {
    
    print("f.lik.con:  begin")
    cat("Initial parameter values:", theta, "\n")
    
  } 
  
  
  if (any(is.na(theta))) {
    
    warning("Problem in f.lik.con: NAs in theta")
    return(NA)
    
  }
  
  if ((length(fct3) > 1) & (length(fct3) != length(x))) {
    
    stop("fct3 incorrect")
    
  }
  
  variance <- 0
  
  for (jj in 1:max(fct3)) variance <- variance + theta[jj] * (fct3 == jj)
  
  regr.par <- theta[(max(fct3) + 1):length(theta)]
  
  
  if (any(!is.finite(theta))) {
    
    theta[!is.finite(theta)] <- par.tmp[!is.finite(theta)]
    
  }
  
  if (sum(theta <= lb) > 0) {
    theta[theta < lb] <- 1.1 * par.tmp[theta < lb]
  }
  
  if (sum(theta >= ub) > 0) {
    theta[theta > ub] <- 0.9 * par.tmp[theta > ub]
  }
  
  expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
      fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, CES = CES, 
      twice = twice, ttt = ttt, y = y, increase = increase, 
      x.mn = x.mn, track = track)
  
  if (any(is.na(expect))) {
    
    warning("NAs in predicted response at parameter values", signif(theta, 4))
    
  }
  
  
  if (dtype %in% c(1, 5, 25, 26)) { 
    # Continuous response 
    
    if(dtype %in% c(1, 5)){
      
      expect <- log(expect)
      yTransformed <- 0
      yTransformed[y > 0] <- log(y)
      
      VdetlimTransformed <- log(Vdetlim)
      cens.upTransformed <- log(cens.up)
      
    } else if (dtype == 26) {
      
      expect <- sqrt(expect)
      yTransformed <- sqrt(y)
      
      VdetlimTransformed <- sqrt(Vdetlim)
      cens.upTransformed <- sqrt(cens.up)
      
    } else {
      
      yTransformed <-  y
      
      VdetlimTransformed <- Vdetlim
      cens.upTransformed <- cens.up
      
    }
    
    score1 <- (y > 0) * (-0.5 * log(2 * pi * variance) - 
          ((yTransformed - expect)^2)/(2 * variance))
    score.detlim <- log(pnorm((VdetlimTransformed - expect)/sqrt(variance)))
    score.detlim[!is.finite(score.detlim)] <- 0
    score.censup <- log(1 - pnorm((cens.upTransformed - expect)/sqrt(variance)))
    score.censup[!is.finite(score.censup)] <- 0
    score2 <- (y == -1000) * score.detlim + (y == -2000) * 
        score.censup
    score <- score1 + score2
    
    
  } else if (dtype %in% c(10, 15, 250, 260)) {
    # Continuous summary response
    
    if(dtype %in% c(10, 15)){
      
      expect <- log(expect)
      
    } else if (dtype == 260){
      
      expect <- sqrt(expect)
      
    }
    
    if (model.ans != 11) {
      
      dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
      
    } else {
      
      dum <- (nn - 1) * sd2.log
      
    }
    
    score <- -(nn * log(sqrt(2 * pi * variance)) + dum/(2 * variance))
    
  } 
  
  if(track)	print("f.lik.con END")
  
  return(-sum(score))
  
  
}



#' Calculate expected response values, for continuous response
#' 
#' This function is used also for quantal response (f.expect.cat)
#' @param model.ans integer, determines the type of model to be fitted
#' @param x numeric vector, the dose values
#' @param regr.par numeric vector, regression parameter values
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param fct3 numeric, value for parameter var
#' @param fct4 numeric, value for parameter c
#' @param fct5 numeric, value for parameter d
#' @param CES numeric, value for the CES
#' @param twice boolean, if some parameter values are equal, see f.execute()
#' @param ttt numeric, time variable
#' @param y numeric vector, the response values
#' @param increase boolean, whether the response values are increasing or 
#' decreasing for increasing dose values 
#' @param x.mn numeric value, the mean of dose values, see f.execute()
#' @param par.start numeric vector, regression parameter values (e.g. MLEs)
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, the expected response values under the estimated 
#' regression model
#' @export
f.expect.con <- function (model.ans, x, regr.par = 0, fct1 = 1, fct2 = 1, fct3 = 1, 
    fct4 = 1, fct5 = 1, CES = NA, twice, ttt = 0, y = 0, 
    increase, x.mn = NA, par.start = NA, track = FALSE) {
  
  if (track)	print("f.expect.con")
  
  nr.aa <- max(fct1)
  nr.bb <- max(fct2)
  nr.var <- max(fct3)
  nr.cc <- max(fct4)
  nr.dd <- max(fct5)
  
  if (model.ans != 11) {
    
    nrp <- length(regr.par)
    aa0 <- rep(0, length(x))
    aa.tmp <- regr.par[1:nr.aa]
    for (ii in (1:nr.aa))
      aa0 <- aa0 + aa.tmp[ii] * (fct1 == ii)
    bb0 <- rep(0, length(x))
    bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
    for (jj in (1:nr.bb))
      bb0 <- bb0 + bb.tmp[jj] * (fct2 == jj)
    par3 <- regr.par[nr.aa + nr.bb + 1]
    if (length(par3) == 0 || is.na(par3)) 
      par3 <- 0
    par4 <- regr.par[nr.aa + nr.bb + nr.cc + 1]
    cc0 <- par3
    dd0 <- par4
    
    if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 19, 
        20, 24, 25, 41, 42, 46)) {
      if (max(fct4) > 1) {
        cc0 <- rep(0, length(x))
        cc.tmp <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + nr.bb + nr.cc)]
        for (kk in (1:nr.cc))
          cc0 <- cc0 + cc.tmp[kk] * (fct4 == kk)
      }
    }
    
    if (model.ans %in% c(3, 8, 13, 18, 23)) {
      dd0 <- par3
      if (max(fct5) > 1) {
        dd0 <- rep(0, length(x))
        dd.tmp <- regr.par[(nr.aa + nr.bb + 1):length(regr.par)]
        for (kk in (1:nr.dd))
          dd0 <- dd0 + dd.tmp[kk] * (fct5 == kk)
      }
    }
    
    if (model.ans %in% c(5, 6, 10, 15, 20, 25, 41, 42, 46)) {
      
      dd0 <- par4
      if (max(fct5) > 1) {
        dd0 <- rep(0, length(x))
        dd.tmp <- regr.par[(nr.aa + nr.bb + nr.cc + 1):length(regr.par)]
        for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
              (fct5 == kk)
      }
      
    }
    
  }
  
  switch(as.character(model.ans), 
      
      # Null model
      "1" = y.expect <- aa0, 
      
      '2' = y.expect <- aa0 * exp(bb0 * x + cc0 * ttt),
      
      '3' = y.expect <- aa0 * exp(bb0 * (x^dd0)),
      
      '4' = y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * x)),
      
      '5' = y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * (x^dd0))),
      
      "11" = { # Full model
        
        x.gr <- levels(factor(x))
        x.fact <- factor(x)
        y.tmp <- rep(NA, length(x))
        y.expect <- rep(0, length(x))
        if (nr.var > 1 & nr.aa == 1 & nr.bb == 1) {
          for (jj in (1:nr.var)) for (ii in (1:length(x.gr))) {
              y.tmp <- y[x.fact == x.gr[ii] & fct3 == jj]
              if (length(y.tmp) > 0) {
                y.mn <- exp(mean(log(y.tmp)))
                y.expect <- y.expect + y.mn * (x.fact == 
                      x.gr[ii]) * (fct3 == jj)
              }
              y.mn <- y.tmp
            }
        } else if (twice) 
          for (jj in (1:nr.bb)) 
            for (ii in (1:length(x.gr))) {
              y.tmp <- y[x.fact == x.gr[ii] & fct2 == jj]
              if (length(y.tmp) > 0) {
                y.mn <- exp(mean(logb(y.tmp)))
                y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * (fct2 == jj)
              }
            } else 
        if (!twice) 
          for (jj in (1:nr.aa))
            for (kk in (1:nr.bb))
              for (ii in (1:length(x.gr))) {
                y.tmp <- y[x.fact == x.gr[ii] & fct1 == jj & fct2 == kk]
                if (length(y.tmp) > 0) {
                  y.mn <- exp(mean(logb(y.tmp)))
                  y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * (fct1 == jj) * (fct2 == kk)
                }
              }
      }, 
      
      "13" = { #E3 - CED
        y.expect <- aa0 * (CES + 1)^((x/bb0)^dd0)
        
      },
      
      '14' = {
        y.expect <- aa0 * (cc0 - (cc0 - 1) * ((CES + 1 - cc0)/(1 - cc0))^(x/bb0))
      },
      
      "15" = { #E5 - CED
        y.expect <- aa0 * (cc0 - (cc0 - 1) * 
              ((CES + 1 - cc0)/(1 - cc0))^((x/bb0)^(dd0)))
        
      }, 
      
      '16' = {
        
        # created in f.start.con if doesn't exist
        CES.16 <- list(CES2 = 0.1, CES1 = 0.05)
        CES2 <- CES.16$CES2
        CES1 <- CES.16$CES1
#				if (par4 <= 1) 
#					print(c("f.expect.con, value of BMD-ratio:", par4))
        dd <- f.uniroot.BMDratio(ratio = par4, cc = cc0, 
            CES1 = CES1, CES2 = CES2, track = track)
        y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd))
        if (is.na(dd)) y.expect <- rep(0, length(x))
        
      },
      
      '17' =  y.expect <- aa0 * (1 - x/(bb0 + x)),
      
      '18' = {
        if (increase == 1) 
          y.expect <- aa0 * (1 - x^dd0/(-(-bb0)^dd0 + x^dd0))
        if (increase == -1)
          y.expect <- aa0 * (1 - x^dd0/((bb0)^dd0 + x^dd0))
      },
      
      '19' = {
        y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(bb0 + x))
      },
      
      '20' = {
        y.expect <- aa0 * (1 + ((cc0 - 1) * x^dd0)/(sign(bb0) * (abs(bb0)^dd0 + x^dd0)))
      },
      
      '21' = {
        
        b2 <- regr.par[nr.aa + nr.bb + 1]
        c1 <- regr.par[nr.aa + nr.bb + 2]
        c2 <- regr.par[nr.aa + nr.bb + 3]
        dd <- regr.par[nr.aa + nr.bb + 4]
        y.expect.1 <- c1 - (c1 - 1) * exp(-(x/bb0)^dd)
        y.expect.2 <- c2 - (c2 - 1) * exp(-(x/b2)^dd)
        y.expect <- aa0 * (y.expect.1 * y.expect.2)
        
      },
      
      "23" = { #H3 - CED
        
        dum <- f.bb.con(model.ans, cc = NA, dd = dd0, CED = bb0, CES, track = track)
        
        if (increase == 1){
          
          y.expect <- aa0 * (1 - x^dd0/(-(-dum)^dd0 + x^dd0))
          
        } else {
          
          y.expect <- aa0 * (1 - x^dd0/((dum)^dd0 + x^dd0))
          
        }
      }, 
      
      '24' = {
        
        dum <- f.bb.con(model.ans, cc = cc0, dd = NA, CED = bb0, CES, track = track)
        y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(dum + x))
        
      },
      
      "25" = { #H5 - CED
        
        dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, CES, track = track)
        y.expect <- aa0 * (1 + ((cc0 - 1) * x^dd0)/(dum^dd0 + x^dd0))
        
      },
      
      '26' = {
        
        y.expect <- cc0 + aa0 * x^bb0
        
      })
  
  if (track)	print("f.expect.con:  END")
  
  return(y.expect)
  
}

#' Determine value for parameter b for continuous response
#' @param model.ans integer, indicates the model that will be fitted;
#' one of \code{c(13, 15, 23, 25)}
#' @param cc numeric, value for parameter c
#' @param dd numeric, value for parameter d
#' @param CED numeric, value for the CED
#' @param CES numeric, value for the CES
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric, value for parameter b
#' @export
f.bb.con <- function (model.ans, cc, dd, CED, CES, track = FALSE) {
  
  if (track)	print("f.bb.con")
  
  switch(as.character(model.ans),
      
      "13" = {
        
        bb <- log(CES + 1)/CED^dd
        
      }, 
      
      "15" = {
        
        CES.tmp <- -abs(CES) * (cc < 1) + CES * (cc >= 1)
        dum <- (CES.tmp + 1 - cc)/(1 - cc)
        bb <- -log(dum)/(CED^dd)
        
      }, 
      
      "23" = {
        
        if (CES >= 0) {
          
          bb <- -CED/(CES/(1 + CES))^(1/dd)
          
        } else {
          
          bb <- CED/(-(CES/(1 + CES)))^(1/dd)
          
        }
        
      }, 
      
      '24' = {
        
        bb <- CED/(CES/(cc - 1 - CES))
        
      },
      
      "25" = {
        
        if (cc == 0) {
          
          bb <- -CED * ((1 + CES)/CES)^(1/dd)
          
        } else {
          
          CES <- abs(CES)
          if (cc >= 1) {
            
            bb <- CED/((CES/(cc - 1 - CES))^(1/dd))
            
          } else {
            
            bb <- CED/((-CES/(cc - 1 + CES))^(1/dd))
            
          }
          
        }
        
      })
  
  if (track)	print("f.bb.con:END")
  
  return(bb)
  
}

#' Determine value for the CED; continuous response
#' @param model.ans integer, indicates the model that will be fitted; 
#' one of \code{c(1, 11, 18, 20, 13, 15, 23, 25)}
#' @param params numeric vector, parameter values for the model parameters a,
#' b, c and d
#' @param CES value for the CES (by default NULL: not used for continuous models)
#' @param dtype response data type, 1 by default
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric, value for the CED 
#' @export
f.inv.con <- function (model.ans, params, CES, dtype = 1, track = FALSE) {
  
  if (track)	print("f.inv.con")
  
  aa <- params[1]
  bb <- params[2]
  par3 <- params[3]
  par4 <- params[4]
  
  switch(as.character(model.ans),
      
      '1' = {
        
        CED <- NA
        warning("For null model CED is not defined")
        
      },
      
      '2' = {
        
        CES.tmp <- if (bb < 0)	-abs(CES)	else	CES
        CED <- (1/bb) * logb(CES.tmp + 1)
        
      },
      
      '3' = {
        
        CES.tmp <- if (bb < 0)	-abs(CES)	else	CES
        CED <- (1/bb) * logb(CES.tmp + 1)
        CED <- CED^(1/par3)
        
      },
      
      '4' = {
        
        CES.tmp <- if (par3 < 1) -abs(CES)	else	CES
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) CED <- -(1/bb) * logb(dum) else {
          if (sum(is.na(dum)) > 0) {
            if(track) warning("parameter c does not allow chosen value for CES")
            CED <- NA
          } else CED <- -(1/bb) * logb(dum)
        }
      },
      
      '5' = {
        
        CES.tmp <- if (par3 < 1) -abs(CES)	else	CES
        dum <- (CES.tmp + 1 - par3)/(1 - par3)
        dum[dum < 0] <- NA
        if (dtype == 3) CED <- (-(1/bb) * logb(dum))^(1/par4) else {
          if (sum(is.na(dum)) > 0) {
            if(track) warning("parameter c does not allow chosen value for CES")
            CED <- NA
          } else CED <- (-(1/bb) * logb(dum))^(1/par4)
        }
        
      },
      
      '11' = CED <- NA,
      
      '13' = CED <- bb,
      
      '15' = CED <- bb,
      
      '17' = {
        
        CES.tmp <- if (bb > 0)	-abs(CES)	else	CES
        CED <- -bb * (CES.tmp/(CES.tmp + 1))
      }, 
      
      '18' = {
        CES.tmp <- if (bb > 0)	-abs(CES)	else	CES
        CED <- ((-sign(bb) * CES.tmp)/(1 + CES.tmp))^(1/par3) * (abs(bb))
      }, 
      
      '19' = {
        CES.tmp <- if (par3 < 1) -abs(CES)	else	CES
        CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))
      }, 
      
      '20' = {
        CES.tmp <- if (par3 < 1) -abs(CES)	else	CES
        if (par3 == 0) 
          CED <- ((-(bb^par4) * CES.tmp)/(1 + CES.tmp))^(1/par4) else 
          CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))^(1/par4)
      },
      
      '23' = CED <- bb,
      
      '25' = CED <- bb
  
  )
  
  
#  } 
#  else if(model.ans == 18){
#    
#    CES.tmp <- CES * (par3 >= 1) - abs(CES) * (par3 < 1)
#    CED <- ((-sign(bb) * CES.tmp)/(1 + CES.tmp))^(1/par3) * (abs(bb))
#    
#  } else if(model.ans == 20){
#    
#    CES.tmp <- CES * (par3 >= 1) - abs(CES) * (par3 < 1)
#    
#    if (par3 == 0){
#      
#      CED <- ((-(bb^par4) * CES.tmp)/(1 + CES.tmp))^(1/par4)
#      
#    } else {
#      
#      CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))^(1/par4)
#      
#    }
#    
  
  # Gives errors, needed?
  if (is.na(CED) && bb < 1e-10) 
    CED <- 1e+10
  
  if (track)	print("f.inv.con : END")
  
  return(CED)
  
}

#' Determine values for the CED; continuous response
#' @param model.ans integer, indicates the model that will be fitted; 
#' one of \code{c(1, 11, 18, 20, 13, 15, 23, 25)}
#' @param regr.par.matr numeric matrix, each row contains a vector of parameter
#' values as needed by the function f.inv.con()
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, calculated CED values
#' @export
f.ced.con <- function (model.ans, regr.par.matr, track = FALSE) {
  
  if (track) 
    print("f.ced.con")
  
  CED <- apply(regr.par.matr, 1, function(par.tmp){
        
        f.inv.con(model.ans, par.tmp, track = track)
        
      })
  
  CED.lst <- list(CED = CED)
  
  if (track) 
    print("f.ced.con: END")
  
  return(CED.lst)
}





#' Fit linear model for continuous response
#' @param xtmp numeric vector, values for the independent variable
#' @param ytmp numeric vector, values for the dependent variable
#' @param model.ans integer, indicates which model will be fitted
#' @param dtype integer, indicates the type of response that is used
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return data frame with intercept and slope of the fitted model, the variable
#' increase is 1 if the slope is positive and -1 otherwise 
#' @export
f.start.lm.con <- function (xtmp, ytmp, model.ans, dtype, track = FALSE) {
  
  if (track)	print("f.start.lm.con")
  
  if (dtype %in% c(25, 250)) {
    
    yTransformed <- ytmp
    
  } else if (dtype %in% c(26, 260)) {
    
    yTransformed <- sqrt(ytmp)
    
  } else {
    
    yTransformed <- log(ytmp)
    
  }
  
  res.lm <- lm(yTransformed ~ xtmp)
  
  intercept <- res.lm[[1]][1]
  slope <- res.lm[[1]][2]
  
  if (dtype %in% c(1, 5, 10, 15)) {
    
    intercept <- exp(intercept)
    
  } else if (dtype %in% c(26, 260)) {
    
    intercept <- intercept^2
    
  } 
  
  increase <- 1 * (slope >= 0) - 1 * (slope < 0)
  
  if (!(model.ans %in% 2:15)) {
    
    yTransformed <- 1/ytmp
    res.lm <- lm(yTransformed ~ xtmp)
    intercept <- 1/res.lm[[1]][1]
    slope <- 1/res.lm[[1]][2]/intercept
    
  }
  
  if (track)	print("f.start.lm.con:  END")
  
  toReturn <- data.frame(intercept = intercept, slope = slope, increase = increase)
  
  return(toReturn)
  
}

#' Calculate residuals: difference of observed and expected response value
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.resid.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.resid.con")
  
  with(ans.all, {
        
        pred.value <- f.expect.con(model.ans, x, regr.par = regr.par, 
            fct1 = fct1, fct2 = fct2, fct3 = fct3, fct5 = fct5, 
            CES = CES, twice = twice, ttt = 0, y = yy, 
            increase = increase, x.mn = x.mn, track = track)
        
        if (dtype %in% c(1, 5, 10, 15)) {
          
          regr.resid <- log(yy) - log(pred.value)
          
        } else if (dtype %in% c(25, 250)) { 
          
          regr.resid <- yy - pred.value 
          
        } else if (dtype %in% c(26, 260)){
          
          regr.resid <- sqrt(yy) - sqrt(pred.value)
          
        } 
        
        ans.all$regr.resid <- regr.resid
        ans.all$pred.value <- pred.value
        
        if (track) 
          print("f.resid.con:  END")
        
        return(ans.all)
        
      })
}



#' Fit the model for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm4.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.mm4.con")
  
  ans.all <- with(ans.all, {
        
        ans.all$conf.int <- matrix(NA, ncol = 2)
        
        if(!is.null(ans.all$parameterConstraints))
          ans.all <- f.constr.con(ans.all, track = track)
        
        if (track) 
          cat("\n\n optimizing .... \n\n")
        
        ans.all <- f.nlminb(ans.all, track = track)
        f.hit.constr(ans.all, track = track)
        
        MLE <- ans.all$MLE
        
        ans.all$regr.par <- MLE[-(1:nr.var)]
        
        if (model.ans != 11){
          
          ans.all$regr.par.matr <- f.pars(ans.all, track = track)$regr.par.matr
          
        } 
        
        ans.all <- f.resid.con(ans.all, track = track)
        
        if (dtype %in% c(5, 15)) {
          if (dtype == 15 || nest.no != 0) {
            if (max(Vdetlim, na.rm = T) != 0) 
              warning("Litter effects not implemented for nonzero detection limit") 
            else if (nr.var > 1) 
              warning("Litter effects not implemented for different within group variances") 
            else ans.all <- f.nested.con(ans.all, track = track)
          } else dtype <- 1
        }
        
        ans.all$boot <- FALSE
        MLE <- ans.all$MLE
        
        if (model.ans %in% c(14, 15, 24, 25))
          cc.OK <- f.check.cc(ans.all, track = track)
        
        covar.ans <- 1
        fit.res <- ans.all$fit.res
        varcov.matr <- NA
        corr.matr <- NA
        if (ans.all$converged) 
          if (length(MLE[MLE == 0]) != 0) {
            warning("The model has too many parameters, use a simpler model")
            covar.ans <- 1
          }
        if (covar.ans == 2) {
          fit.res$lower <- ans.all$lb
          fit.res$upper <- ans.all$ub
          fit.res$scale <- scale.dum
          varcov.matr <- vcov(fit.res)
          cat("\n\n variance-covariance matrix: \n")
          print(varcov.matr)
          if (varcov.matr[1] == -1000) {
            cat("\ntry fitting again, with scale adjusted\n")
            corr.matr <- 1000
          }
          else {
            tmp.matr <- sqrt(diag(varcov.matr))
            tmp.matr <- diag(1/tmp.matr)
            corr.matr <- tmp.matr %*% varcov.matr %*% tmp.matr
            cat("\n correlation matrix: \n")
            print(corr.matr)
            cat("\n maximum correlation:     ", round(max(abs(lower.tri(corr.matr) * 
                                corr.matr)), 4), "\n")
          }
        }
        
        ans.all$par.start <- ans.all$MLE
        ans.all$list.logic <- FALSE
#        if (model.ans %in% c(12:15, 22:25)) {
#          
#          ans.all$CED <- f.ced.con(model.ans, ans.all$regr.par.matr, track = track)
#          
#        }
        
        if (track) 
          print("f.mm4.con: END")
        
        return(ans.all)
        
      })        
  
  
  
}


#' Calculate the CED values for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm6.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.mm6.con")
  
  ans.all <- with(ans.all, {
        
        if (fitted || list.logic) {
          
          ans.all <- f.pars(ans.all, track = track)
          
        }
        
        CED <- f.ced.con(model.ans, regr.par.matr, track = track)$CED
        
        CED.unscaled <- CED
        
        if (max(nr.aa, nr.bb) == 1) {
          
          CED.unscaled <- sf.x * CED[1]
          
        } else {
          
          CED.unscaled <- sf.x * CED
          
        }
        
        
        if (dtype %in% c(5, 15)) {
          
          ans.all$CED <- CED
          return(ans.all)
          
        }
        
        f.check.cc(ans.all, track = track)
        
        ans.all$group <- 0
        ans.all <- f.CI(ans.all, track = track)
        
        
        nr.CED <- max(fct2)
        
        if (!all(is.na(conf.int))) {
          
          CED.unique <- unique(CED)
          
          if (length(CED.unique) > 1) {
            
            CED.unscaled <- sf.x * CED.unique[nr.CED]
            
          } else {
            
            CED.unscaled <- sf.x * CED.unique
            
          }          
        }
        
        
        ans.all$CED <- CED
        ans.all$CED.origScale <- CED.unscaled
        ans.all$BMR <- NA
        
        if (track)	print("f.mm6.con:  END")
        
        return(ans.all)
        
      })
}

#' Define lower and upper bounds for the model parameters; continuous response 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.constr.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.constr.con")
  
  with(ans.all, {        
        
        # MV added this
        if(!is.null(ans.all$parameterConstraints)){
          
          lower <- ans.all$parameterConstraints[,"lowerBound"]
          upper <- ans.all$parameterConstraints[,"upperBound"]
          
          lower.var <- lower[1]
          lower.aa <- lower[2]
          lower.bb <- lower[3]
          upper.var <- upper[1]
          upper.aa <- upper[2]
          upper.bb <- upper[3]
          
          if (model.ans %in% c(15, 25)) {
            
            lower.cc <- lower[4]
            upper.cc <- upper[4]
            lower.dd <- lower[5]
            upper.dd <- upper[5]
            
          } else if (model.ans %in% c(13, 23)) {
            
            lower.dd <- lower[4]
            upper.dd <- upper[4]
            
          }
          
        } else {
          
          lower.var <- 1e-06
          upper.var <- 10
          
          if (dtype %in% c(25, 250)) {
            
            upper.var <- 1e+10
            
          } else if (dtype %in% c(26, 260)) {
            
            upper.var <- 1e+04
            
          }           
          
          if (model.ans %in% c(1:10, 12:20, 22:25, 47)) {
            
            nr.aa <- max(fct1)
            nr.bb <- max(fct2)
            nr.cc <- max(fct4)
            nr.dd <- max(fct5)
            nr.var <- max(fct3)
            
            if (cont) {
              
              lower.aa <- min(yy / 100)
              upper.aa <- max(yy * 100)
              
            } else 
            if (model.type == 2) {
              
              lower.aa <- 1e-06
              upper.aa <- 10
              
            }
            
            lower.bb <- -Inf
            upper.bb <- Inf
            
            if (model.ans %in% c(4:6, 9:10, 13:15, 16, 19:20, 23:25, 47)) {
              
              lower.bb <- 1e-06
              
            }
            
            if (model.ans == 3) {
              
              if (increase == -1) {
                lower.bb <- -Inf
                upper.bb <- 0
              }
              if (increase == 1) {
                lower.bb <- 0
                upper.bb <- Inf
              }
              
            }
            
            if (increase == 1) {
              
              lower.cc <- 1.1
              upper.cc <- 1e+06
              
            } else {
              
              lower.cc <- 1e-06
              upper.cc <- 0.9
              
            }
            
            lower.dd <- 0.01
            upper.dd <- 10
            
            if (model.ans == 1) {
              
              lower <- lower.aa
              upper <- upper.aa
              
            } else if (model.ans %in% c(2, 7, 12, 17, 22)) {
              
              lower <- c(lower.aa, lower.bb)
              upper <- c(upper.aa, upper.bb)
              
            } else if (model.ans %in% c(3, 8, 13, 18, 23)) {
              
              lower <- c(lower.aa, lower.bb, lower.dd)
              upper <- c(upper.aa, upper.bb, upper.dd)
              
            } else if (model.ans %in% c(4, 9, 14, 19, 24)) {
              
              lower <- c(lower.aa, lower.bb, lower.cc)
              upper <- c(upper.aa, upper.bb, upper.cc)
              
            } else if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25)) {
              
              lower <- c(lower.aa, lower.bb, lower.cc, lower.dd)
              upper <- c(upper.aa, upper.bb, upper.cc, upper.dd)
              
            } else if (model.ans == 47) {
              
              lower.qq <- 1e-12
              upper.qq <- 100
              lower <- c(lower.aa, lower.bb, lower.qq, lower.dd)
              upper <- c(upper.aa, upper.bb, upper.qq, upper.dd)
              
            }
            
            if (cont) {
              
              lower <- c(lower.var, lower)
              upper <- c(upper.var, upper)
              
            } else {
              
              lower <- c(NA, lower)
              upper <- c(NA, upper)
              
            }
            
          } else if (model.ans == 11) {
            
            lower <- c(lower.var, regr.par)
            upper <- c(upper.var, regr.par)
            
          }
          
        }        
        
        lower.var <- lower[1]
        lower.aa <- lower[2]
        lower.bb <- lower[3]
        
        if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 24, 25, 26)) {
          
          lower.cc <- lower[4]
          upper.cc <- upper[4]
          lower.dd <- lower[5]
          upper.dd <- upper[5]
          
        } else if (model.ans %in% c(3, 8, 13, 18, 23)) {
          
          lower.dd <- lower[4]
          upper.dd <- upper[4]
          
        }
        
        upper.var <- upper[1]
        upper.aa <- upper[2]
        upper.bb <- upper[3]
        
        
        # Repeat lower and upper bounds nr.<par> times
        if (model.ans == 1) {
          lb <- c(rep(lower.aa, nr.aa))
          ub <- c(rep(upper.aa, nr.aa))
        } else if (model.ans %in% c(2, 7, 12, 17, 22, 35, 40, 44)) {
          lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb))
          ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb))
        } else if (model.ans %in% c(3, 8, 13, 18, 23)) {
          lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.dd, nr.dd))
          ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.dd, nr.dd))
        } else if (model.ans %in% c(4, 9, 14, 19, 24, 26, 27, 30)) {
          lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.cc, nr.cc))
          ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.cc, nr.cc))
        } else if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 28, 29, 
            31, 34, 36, 37, 41, 42, 43)) {
          lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.cc, nr.cc), rep(lower.dd, nr.dd))
          ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.cc, nr.cc), rep(upper.dd, nr.dd))
        } else if (model.ans == 11) {
          lb <- regr.par
          ub <- regr.par
        } else if (model.ans == 21) {
          lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd, lower.ee, rep(lower.ff, nr.dd))
          ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd, upper.ee, rep(upper.ff, nr.dd))
        }
        
        if (cont) {
          
          lb <- c(rep(lower.var, nr.var), lb)
          ub <- c(rep(upper.var, nr.var), ub)
          
        } else {
          
          par.start <- c(regr.start, th.start, sig.start)
          nrp <- length(regr.start)
          npar <- length(par.start)
          if (model.ans %in% c(2:3, 7:8)) 
            ub[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
          if (model.ans %in% c(4:5, 9:10)) 
            lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
          if (model.ans %in% c(17:20)) 
            lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
          if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 24, 25)) {
            lb[nr.aa + nr.bb + 1] <- 1e-06
            ub[nr.aa + nr.bb + 1] <- 0.999999
          }
          if (model.ans %in% c(12:15, 22:25)) {
            if (nr.aa > 1 && nr.bb == 1) {
              lb[2:(nr.aa + 1)] <- eps
              ub[2:(nr.aa + 1)] <- Inf
            }
            else {
              lb[(nr.aa + 1):(nr.aa + nr.bb)] <- eps
              ub[(nr.aa + 1):(nr.aa + nr.bb)] <- Inf
            }
          }
          lb[nrp + 1] <- th.0.start
          ub[nrp + 1] <- th.0.start
          if (nth > 1) {
            lb[(nrp + 2):npar] <- -Inf
            ub[(nrp + 2):npar] <- 1e-10
          }
          lb[npar] <- sig.start
          ub[npar] <- sig.start
          if (max(as.numeric(fct3)) > 1) {
            lb[npar] <- 0
            ub[npar] <- Inf
          }
          if (model.ans %in% c(3, 5, 18, 20, 13, 15, 23, 25)) 
            if (constr != -Inf) 
              lb[length(par.start) - nth - 1] <- constr
          if (model.ans %in% c(4, 5, 9, 10, 14, 15, 24, 25)) {
            lb[nr.aa + nr.bb + 1] <- 0
            ub[nr.aa + nr.bb + 1] <- 1
          }
          if (dtype == 6) {
            text.par <- c("alfa", text.par)
            par.start <- c(10, par.start)
            if (dtype == 6) 
              if (!(model.ans == 14 & model.type == 1)) 
                par.start[1] <- ans.all$alfa.start
            npar <- length(par.start)
            lb <- c(1e-10, lb)
            ub <- c(Inf, ub)
          }
        }
        
        ans.all$lb <- lb
        ans.all$ub <- ub
        ans.all$lower <- lower
        ans.all$upper <- upper
        
        if (track) 
          print("f.constr.con:  END")
        
        return(ans.all)
        
      })
}


#' compute the root of the BMD ratio?
#' 
#' Used in f.expect.con for model.ans == 24 (quantal model)
#' @param ratio ratio of interest
#' @param cc value for c parameter
#' @param CES1 first value for CES
#' @param CES2 second value for CES
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return root
#' @export
f.uniroot.BMDratio <- function (ratio, cc, CES1, CES2, track = FALSE) 
{
  if (track)	print("f.uniroot.BMDratio")
  
  if (ratio <= 1) {
    cat("ATTENTION (f.uniroot):  BMD ratio is smaller than (or equal to) one\n")
    return(NA)
  }
  
  if ((cc > 1 - CES2) & (cc < 1 + CES2)) {
    cat("ATTENTION: parameter c (', cc, ') does not allow CES2 of:", CES2, "\n")
    if (cc > 1)	cc <- 1.0000001 * (1 + CES2)
    if (cc < 1)	cc <- 0.9999999 * (1 - CES2)
  }
  
  if (cc < 1) {
    CES1 <- -CES1
    CES2 <- -CES2
  }
  
  f.c <- function(dd, ratio, cc, CES1, CES2) {
    yy <- (log((CES2 + 1 - cc)/(1 - cc))/log((CES1 + 1 - cc)/(1 - cc)))^(1/dd)
    return(yy - ratio)
  }
  
  dd.low <- 0.1
  dd.upp <- 100
  try1 <- f.c(dd.low, ratio, cc, CES1, CES2)
  try2 <- f.c(dd.upp, ratio, cc, CES1, CES2)
  
  if (is.na(try1)) 
    return(NA)
  if (is.na(try2)) 
    return(NA)
  if (sign(try1) == sign(try2)) {
    cat("\nATTENTION:  no root found for parameter d")
    return(NA)
  }
  
  root.out <- uniroot(f.c, interval = c(dd.low, dd.upp), ratio = ratio, 
      cc = cc, CES1 = CES1, CES2 = CES2)
  
  return(root.out$root)
  
}


#' Calculations for nested factors (clustered response data)
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.nested.con <- function (ans.all, track = FALSE) {
  
  if (track) print("f.nested.con")
  
  with(ans.all, {
        
        if (dtype %in% c(5, 15)) {
          
          if (dtype == 5) {
            fact.nest <- data.0[, nest.no]
            lev <- factor(fact.nest)
            nest.size <- as.numeric(tapply(fact.nest, fact.nest, length))
            no.of.nests <- length(as.numeric(nest.size))
          }
          
          if (dtype == 15) {
            nest.size <- nn
            no.of.nests <- length(as.numeric(nest.size))
            lev <- factor(1:no.of.nests)
          }
          
          aov.res <- aov(regr.resid ~ lev)
          aov.summ <- summary(aov.res)
          
          # MV added condition because results in error for dtyp = 15
          if(dtype == 5)
            ans.all$p.value <- aov.summ[[1]][[5]][1]
          
          MS <- aov.summ[[1]][[3]]
          df <- aov.summ[[1]][[1]]
          
          MS.between <- MS[1]
          MS.within <- MS[2]
          df1 <- df[1]
          df2 <- df[2]
          Q1 <- 0
          Q2 <- 0
          Q3 <- 0
          Q4 <- 0
          dum.lev <- as.numeric((lev))
          nij <- as.numeric(tapply(dum.lev, dum.lev, length))
          Q1 <- sum(nij)
          Q2 <- sum(nij^2)
          n0 <- 1/df1 * (Q1 - Q2/Q1)
          inter.var <- var(regr.resid) - MS.within
          
          if (inter.var <= 0) {
            
            warning("ATTENTION: variance between clusters was found to be nonpositive\n              and will be set to zero, i.e., clustering is omitted")
            ans.all$dtype <- 1
            return(ans.all)
            
          }
          
          var.between.2 <- (MS.between - MS.within)/n0
          intra.var <- signif(MS.within, 4)
          
        }
        
        if (dtype == 15) {
          
          no.of.nests <- length(nn)
          SS <- sum(nn * sd2.log)
          df1 <- length(nn) - 1
          df2 <- sum(nn) - length(nn)
          intra.var <- SS/df2
          inter.var <- MLE[1] - intra.var
          
          ans.all$nest.size <- nn
          
        }
        
        if (dtype == 5) {
          
          ans.all$fact.nest <- fact.nest
          ans.all$no.of.nests <- no.of.nests
          ans.all$nest.size <- nest.size
          
        }
        
        ans.all$inter.var <- inter.var
        ans.all$intra.var <- intra.var
        ans.all$df1 <- df1
        ans.all$df2 <- df2
        
        if (track) print("f.nested.con:  END")
        return(ans.all)
        
      })
}



#' Main function for bootstrap method for calculating CI of CED
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm7.con <- function (ans.all, track = track) {
  
  if (track) print("f.mm7.con")
  
  with(ans.all, {
        
        ans.all.bt <- ans.all
        
        y.true <- f.expect.con(model.ans, x, MLE[(max(fct3) + 1):length(MLE)], 
            fct1 = fct1, fct2 = fct2, fct3 = fct3, 
            fct5 = fct5, CES = CES, twice = twice, y = yy, increase = increase)
        
        
        ans.all.bt$y.true <- y.true
        ans.all.bt$x <- x
        ans.all.bt$fct1 <- fct1
        ans.all.bt$fct2 <- fct2
        ans.all.bt$fct3 <- fct3
        ans.all.bt$fct5 <- fct5
        
        if (dtype %in% c(10, 15, 110, 250, 260)) {
          
          ans.all.bt$x <- numeric()
          ans.all.bt$fct1 <- numeric()
          ans.all.bt$fct2 <- numeric()
          ans.all.bt$fct3 <- numeric()
          ans.all.bt$fct5 <- numeric()
          ans.all.bt$y.true <- numeric()
          
          for (ii in 1:length(x)) {
            
            ans.all.bt$y.true <- c(ans.all.bt$y.true, rep(y.true[ii], 
                    nn[ii]))
            ans.all.bt$x <- c(ans.all.bt$x, rep(x[ii], nn[ii]))
            ans.all.bt$fct1 <- c(ans.all.bt$fct1, rep(fct1[ii], 
                    nn[ii]))
            ans.all.bt$fct2 <- c(ans.all.bt$fct2, rep(fct2[ii], 
                    nn[ii]))
            ans.all.bt$fct3 <- c(ans.all.bt$fct3, rep(fct3[ii], 
                    nn[ii]))
            ans.all.bt$fct5 <- c(ans.all.bt$fct5, rep(fct5[ii], 
                    nn[ii]))
            
          }
          
          if (dtype == 10) 
            ans.all.bt$dtype <- 1
          if (dtype == 15) 
            ans.all.bt$dtype <- 5
          
        }
        
        sigma <- rep(0, length(ans.all.bt$x))
        
        for (jj in (1:max(ans.all.bt$fct3))) 
          sigma <- sigma + sqrt(MLE[jj]) * (ans.all.bt$fct3 == jj)
        ans.all.bt$sigma <- sigma
        ans.all.bt$df.corr <- length(ans.all.bt$x) - nrp
        ans.all.bt$par.start <- MLE
        CED.boot <- matrix(nrow = nruns, ncol = length(CED))
        par.boot <- matrix(nrow = nruns, ncol = length(MLE))
        converged.boot <- numeric()
        
        if(track)
          cat("\n\nCalculating confidence interval by bootstrapping .......\n")
        
        
        for (run in 1:nruns) {
          
          set.seed(run)
          boot.lst <- f.boot.con(ans.all.bt)
          par.boot[run, ] <- boot.lst$MLE
          CED.boot[run, ] <- boot.lst$CED
          converged.boot[run] <- boot.lst$conv
          if(track)
            cat(run, "   *Parameter estimates:", signif(par.boot[run, ], 3), 
                "\n     *CED estimate: ", signif(sf.x * CED.boot[run, ], 5), "\n")
        }
        
        ans.all$CED <- CED
        ans.all$par.boot <- par.boot
        ans.all$CED.boot <- CED.boot
        ans.all$converged.boot <- converged.boot
        
        nruns.intend <- length(par.boot[, 1])
        nruns.conv <- sum(converged.boot)
        ans.all$nruns.conv <- nruns.conv
        ans.all$nruns.intend <- nruns.intend
        
        pct.lst <- f.ced.pct(CED = CED, CED.boot = CED.boot, nruns.intend = nruns.intend, 
            sf.x = sf.x, conf.lev = ans.all$conf.lev)
        
        ans.all$conf.int <- pct.lst$conf.int
        
        ans.all$boot <- TRUE
        
        if (track) print("f.mm7.con:  END")
        
        
        return(ans.all)
      })
}


#' One bootstrap run for continuous response
#' @param ans.all.bt list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, estimated MLEs, CED and whether the model converged
#' @export
f.boot.con <- function (ans.all.bt, track = FALSE) {
  
  if (track) print("f.boot.con")
  
  with(ans.all.bt, {
        
#        options(warn = -1)
        noise <- numeric()
        
        if (dtype %in% c(5, 15)) {
          
          inter.sd <- sqrt((inter.var * df1)/rchisq(1, df1))
          intra.sd <- sqrt((intra.var * df2)/rchisq(1, df2))
          if (intra.sd < sqrt(intra.var)) 
            intra.sd <- sqrt(intra.var)
          eps <- rnorm(no.of.nests, 0, inter.sd)
          
          for (ii in 1:no.of.nests) {
            noise <- c(noise, rnorm(nest.size[ii], eps[ii], 
                    intra.sd))
            ans.all.bt$dtype <- 1
          }
          
        } else {
          
          sigma <- sigma * sqrt(df.corr/rchisq(1, df.corr))
          noise <- rnorm(length(x), 0, sigma)
          
        }
        
        
        if (dtype == 25 | dtype == 250) 
          y <- y.true + (noise)
        else if (dtype == 26 | dtype == 260) 
          y <- (sqrt(y.true) + (noise))^2
        else y <- y.true * exp(noise)
        
        
        if (dtype == 250) 
          ans.all.bt$dtype <- 25
        if (dtype == 260) 
          ans.all.bt$dtype <- 26
        ans.all.bt$y <- y
        ans.all.bt$yy <- y
        
        ans.all.bt <- f.nlminb(ans.all = ans.all.bt, track = track)
        par.bt <- ans.all.bt$MLE
        regr.par.bt <- par.bt[(max(fct3) + 1):length(par.bt)]
        ans.all.bt <- f.pars(ans.all.bt)
        regr.bt.matr <- ans.all.bt$regr.par.matr
        
        CED.lst <- f.ced.con(model.ans = model.ans, regr.par.matr = regr.bt.matr, track = track)
        CED.bt <- CED.lst$CED
        
        if (track) print("f.boot.con:  END")
        
        
        return(list(MLE = par.bt, CED = CED.bt, conv = ans.all.bt$converged))
        
      })
}



#' Calculate CI based on bootstrap estimates and confidence level
#' @param CED numeric, the estimated value for CED 
#' @param CED.boot numeric vector, the estimated CED values following bootstrap
#' @param nruns.intend integer, the intended number of bootstrap runs
#' @param sf.x numeric, scaling factor for x; default value is NA
#' @param conf.lev numeric, the confidence level for confidence intervals; 
#' default value is 0.9
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, the confidence interval values and the confidence level as 
#' provided in conf.lev
#' @export
f.ced.pct <- function (CED, CED.boot, nruns.intend, sf.x = NA, conf.lev = 0.9,
    track = FALSE) {
  
  if (track) print("f.ced.pct")
  
  if (nruns.intend < 5) {
    
    cat("\n use larger number of bootstrap runs\n")
    ced.lst <- list(conf.int = matrix(NA, ncol = 2), conf.lev = conf.lev)
    return(ced.lst)
    
  }
  
  CED.nr <- length(CED.boot[1, ])
  CED <- signif(CED, 4)
  tb <- "\t"
  Lmed <- numeric()
  Llow <- numeric()
  Lupp <- numeric()
  conf.int <- matrix(NA, nrow = CED.nr, ncol = 2)
  
  for (jj in (1:CED.nr)) {
    
    CED.bt <- CED.boot[, jj]
    CED.bt <- CED.bt[is.finite(CED.bt)]
    x.lab <- "CED animal"
    if (CED.nr > 1) 
      x.lab <- paste(x.lab, "group", jj)
    
#    if (trace.plt) {
#      f.graph.window(2)
#      hist(CED.bt, xlab = x.lab, nclass = 15, main = "")
#      hist.out <- hist(CED.bt, plot = F, nclass = 15)
#      y.max <- max(hist.out$counts)
#      if (is.R()) 
#        mtext(show, 4, 1, adj = 0, padj = 1, las = 1, 
#            at = y.max, cex = 0.8)
#      else mtext(show, 4, 0.3, adj = 0, las = 1, at = y.max, 
#            cex = 0.8)
#      hist(log10(CED.bt), xlab = paste("log10 of", x.lab), 
#          nclass = 15, main = "")
#      hist.out <- hist(log10(CED.bt), plot = F, nclass = 15)
#      y.max <- max(hist.out$counts)
#    }
    
    Lmed[jj] <- quantile(CED.bt, 0.5)
    Llow[jj] <- quantile(CED.bt, (1 - conf.lev)/2)
    Lupp[jj] <- quantile(CED.bt, 1 - (1 - conf.lev)/2)
    conf.int[jj, ] <- c(sf.x * Llow[jj], sf.x * Lupp[jj])
#    CED.txt <- paste("L50", fct2.txt[jj], sep = "-")
#    
#    if (trace.cat) {
#      cat("\n", x.lab, "\n")
#      cat("5th percentile", tb, "50th percentile", tb, 
#          "95th percentile\n")
#      cat(fct2.txt[jj], tb, sf.x * Llow[jj], tb, tb, sf.x * 
#              Lmed[jj], tb, tb, sf.x * Lupp[jj], "\n")
#    }
#    
#    Llow.txt <- paste("L", 100 * ((1 - conf.lev)/2), sep = "")
#    Lupp.txt <- paste("L", 100 * (1 - (1 - conf.lev)/2), 
#        sep = "")
#    show.tmp <- paste(y.leg, "\nCES", round(CES, digits = 2), 
#        "\n", CED.txt, ": ", Lmed[jj], "\n", Llow.txt, ": ", 
#        Llow[jj], "\n", Lupp.txt, ": ", Lupp[jj], "\nruns:", 
#        nruns.intend, "\nconv:", nruns.conv)
#    if (trace.plt) {
#      if (is.R()) 
#        mtext(show.tmp, 4, 1, adj = 0, padj = 1, las = 1, 
#            at = y.max, cex = 0.8)
#      else mtext(show.tmp, 4, 0.3, adj = 0, las = 1, at = y.max, 
#            cex = 0.8)
#    }
  }
  
  ced.lst <- list(conf.int = conf.int, conf.lev = conf.lev)
#  if (trace.plt) {
#    cat("\n\nnumber of intended bootstrap runs: ", nruns.intend)
#    cat("\nnumber of converged bootstrap runs: ", nruns.conv)
#    cat("\ncritical effect size: ", CES)
#  }
  
  if (track) print("f.ced.pct: END")
  
  
  return(ced.lst)
}

