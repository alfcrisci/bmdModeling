#' Main function for calculations with categorical response type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all 
#' @export
f.cat <- function (ans.all, track = FALSE){
  
  if (track)	print("f.cat")
  
#	# Note: contrary to f.cont, f.start.X not included at the beginning of the function
#	# added, otherwise par.start is missing
#	ans.all <- f.start.bin(ans.all, track = track) # user f.start.bin because model.type == 1
#  
  
  if(!is.null(ans.all$startValues))
    ans.all$main.ans <- c(3, ans.all$main.ans)
  
  for(main.ans.single in c(2, ans.all$main.ans, 15)){
    
    switch(as.character(main.ans.single), 
        
        # 2: define initial values
        '2' = {
          
          ans.all <- with(ans.all, {
                
                if (dtype == 6 && !(ans.all$model.type == 1 & ans.all$model.ans == 14)) 
                  ans.all <- f.dtype6.mn(ans.all, track = track)
                if (dtype == 4) {
                  
                  model.ans.0 <- ans.all$model.ans
                  model.type.0 <- ans.all$model.type
                  fct1.0 <- ans.all$fct1
                  fct2.0 <- ans.all$fct2
                  ans.all$model.ans <- 14
                  ans.all$model.type <- 1
                  
                  if(track)
                    cat("\ncalculating group means using full model\n")
                  
                  ans.all <- f.start.bin(ans.all, track = track)
                  
                  ans.all$model.ans <- model.ans.0
                  ans.all$model.type <- model.type.0
                  ans.all$fct1 <- fct1.0
                  ans.all$fct2 <- fct2.0
                  
                }
                
                # From f.choose.model()
#                full LVM model not allowed, model 14 from classical models will be chosen
                if (model.type == 2 & model.ans == 11) {
                  cat("\nfull LVM model not allowed, model 14 from classical models will be chosen\n")
                  ans.all$model.ans <- 14
                  ans.all$model.type <- 1
                  ans.all$modelname <- "Model 14 (E4 - CED)"
                }
                
                if (dtype == 3) {
                  ans.all$CES <- 0
                  ans.all$ces.ans <- 1
                }
                
                dum1 <- ((model.type == 1) & (model.ans %in% c(15:26, 30:31, 33)))
                dum2 <- ((model.type == 2) & (model.ans %in% c(12:15, 22:25)))
                ans.all$CED.model <- dum1 | dum2
                # end from f.choose.model()                
                
                switch(ans.all$model.type, 
                    ans.all <- f.start.bin(ans.all, track = track),
                    ans.all <- f.start.cat(ans.all, track = track)
                )
                
                if (is.na(ans.all$loglik.old) | ans.all$loglik.old <= -2e+10) {
                  ans.all$errorAdjustStartValues <- TRUE
#                    switch(ans.all$model.type, 
#                        ans.all <- f.start.bin(ans.all, adjust = T), 
#                        ans.all <- f.start.cat(ans.all, adjust = T))
                }
                
                ans.all.tmp <- ans.all
                
              }, return(ans.all.tmp))
          
        },
        
        # 3: change start values
        "3" = { 
          
          if (ans.all$model.type == 1) ans.all <- f.start.bin(ans.all, 
                adjust = TRUE, track = track)
          if (ans.all$model.type == 2) ans.all <- f.start.cat(ans.all, 
                adjust = TRUE, track = track)
          
        },
        
        # 4: fit regression
        '4' = {
          
          ans.all <- f.mm4.cat(ans.all, track = track)
          
        },
        
        # 6: Calculate BMD/CED
        '6' = {
          
          if (!ans.all$fitted){
            
            stop("First fit the model")
            
          }
          
          if(ans.all$model.ans %in% c(1, 11)){
            
            stop("BMD not defined for null or full model")
            
          }
          
          ans.all <- f.mm6.cat(ans.all, track = track)
          
        }, 
        
        # end session
        
        '15' = {
          
          return(ans.all)
          
        })
    
  }
  
}

#' Fit the model for a categorical response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm4.cat <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.mm4.cat")
  
  with(ans.all, {
        
        # MV added this
        if(!is.null(ans.all$parameterConstraints)){
          
          lb <- ans.all$parameterConstraints[,"lowerBound"]
          ub <- ans.all$parameterConstraints[,"upperBound"]
          
          ans.all$lb <- lb
          ans.all$ub <- ub
          
        }
        
        ans.all <- f.nlminb(ans.all, track = track)
        f.hit.constr(ans.all)
        ans.all$CED.matr <- matrix(NA, ncol = 2, nrow = 1)
        loglik <- ans.all$loglik
        MLE <- ans.all$MLE
        
        if (model.type == 1 & model.ans == 14) {#model.type == 1 & 
          
          pi.full <- MLE
          if (dtype == 6) {
            alfa.start <- MLE[1:alfa.length]
            pi.full <- pi.full[-(1:alfa.length)]
          }
          
        }
        
        converged <- ans.all$converged
        aa <- MLE[1:max(fct1)]
        if ((min(x) == 0) & (model.type == 1) & (prod(aa) == 0)) {#& (model.type == 1)
          dum <- (y[x == 0] > 0)
          warning("ATTENTION: one of parameters a is estimated to be zero,",
              "which is impossible for nonzero observed response at dose zero.",
              "Therefore you should refit the model with positive lower bound for parameter a")
        }
        
        regr.par <- if (dtype == 6)	MLE[-(1:alfa.length)]	else	MLE
        
        if (nrp == 0)	regr.par <- numeric()
        
        if (model.type == 2) {
          par.lst <- f.split.par(MLE, nrp, nth, dtype)
          regr.par <- par.lst$regr.par
          th.par <- par.lst$th.par
          sig.par <- par.lst$sig.par
        }
        
        if (model.type == 2 & dtype %in% c(4, 6)) {
          pi.tmp <- f.expect.cat(model.type, model.ans, x, 
              regr.par, th.par, sig.par, fct1, fct2, CES.cat = CES.cat, 
              CES = CES, ttt = ttt, twice = twice, fct3 = fct3, 
              ces.ans = ces.ans, decr.zz = decr.zz, dtype = dtype, 
              fct3.ref = fct3.ref, track = track)
          y.dum <- tapply(pi.tmp, x, mean)
          x.dum <- tapply(x, x, mean)
        }
        
        ans.all$pi.full <- pi.full
        ans.all$regr.par <- regr.par
        ans.all$th.par <- th.par
        ans.all$sig.par <- sig.par
        ans.all$l.ty <- 1
        ans.all <- f.pars(ans.all, track = track)
        
        # remove all plots
        
        if (dtype %in% 2:3) {
          ans.all$combi <- FALSE
          ans.all$plot.type <- 1 # original: 5
          ans.all$categ.ans <- 0
        }
        
        if (model.type == 1 && !model.ans %in% 25:28) 
          if (any(MLE[1:nr.aa] < 0)) 
            warning("The parameter a is estimated to be negative \n",
                "Use maximum likelihood with constraints and refit the model")
        
        ans.all$odt <- odt
        ans.all$MLE <- MLE
        ans.all$alfa.start <- alfa.start
        ans.all$par.start <- MLE
        ans.all$categ.ans <- 1
        ans.all$fitted <- TRUE
        
        if (track)	print("f.mm4.cat:  END")
        
        return(ans.all)
        
      })
  
}

#' Calculate the CED values for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm6.cat <- function (ans.all, track = FALSE) 
{
  if (track)	print("f.mm6.cat")
  
  with(ans.all, {
        
        if (twice) 
          max.lev <- max(max(fct1), max(fct2))
        if (!twice) 
          max.lev <- max(fct1) * max(fct2)
        
        conf.int <- matrix(NA)
        
        if (!fitted) {
          
          if(track)	cat("\nATTENTION: you did not yet fit the model!\n\n")
          regr.par <- regr.start
          th.par <- c(0, th.start)
          sig.par <- sig.start
        }
        
        profile.ans <- ifelse(fitted & CED.model, 1, 2)
#		profile.ans <- menu(c("yes", "no"), title = "\nDo you want to calculate confidence intervals for BMD?")
        
        if (dtype == 3) {
          ces.ans <- ans.all$ces.ans <- 1
          ans.all$CES <- 0
        }
        
        # ask for CES value (done in the UI)
#			if (dtype != 3 && !CED.model) {
#				ces.ans <- menu(c("ED50", "Additional risk, i.e. P[BMD] - P[0]", 
#									"Extra risk, i.e. (P[BMD]-P[0])/(1-P[0])", "", 
#									""), title = "\nWhat type of Benchmark response do you want to consider?")
#				ans.all$ces.ans <- ces.ans
#				switch(ces.ans, CES <- 0, CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of additional risk > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of extra risk > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of difference in z-score > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the percent change in risk > ")))
#			}
        
        ans.all$regr.par <- regr.par
        ans.all$th.par <- th.par
#		ans.all$CES <- CES
        CED.list <- f.ced.cat(ans.all, track = track)
        CED.matr <- CED.list$CED.matr
        
        gr.txt <- CED.list$gr.txt
        ans.all$gr.txt <- gr.txt
        ans.all$response.matr <- CED.list$response.matr
        pi.bg.matr <- CED.list$pi.bg.matr
        
        ans.all$CED.matr <- CED.matr
        
        nr.CED <- length(CED.matr[, 1])
        
        if (profile.ans == 1) {
          
          ans.all$trace <- TRUE
          ans.all$trace.plt <- TRUE
          ans.all$group <- 0
          
          if (model.type == 1 && model.ans == 25 && ces.ans == 1) #model.type == 1 && 
            ans.all$group <- 1:nr.aa
          ans.all <- f.CI(ans.all, track = track)
          
          if (ans.all$update) {
            MLE <- ans.all$MLE
            switch(model.type, dum <- 1, dum <- 0)
            CED.list <- f.ced.cat(ans.all, track = track)
            CED.matr <- CED.list$CED.matr
            ans.all$regr.par <- MLE
            if (dtype == 6) 
              ans.all$regr.par <- ans.all$regr.par[-1]
            if (model.type == 2) {
              par.lst <- f.split.par(MLE, nrp, nth, dtype)
              ans.all$regr.par <- par.lst$regr.par
              ans.all$th.par <- c(0, par.lst$th.par)
              ans.all$sig.par <- par.lst$sig.par
            }
          }
          
          ans.all$CED.matr <- CED.matr
          
          conf.int <- ans.all$conf.int
          
          nr.CED <- length(conf.int[, 1])
          
#			if (nr.CED > 1) 
#				for (ii in 1:nr.CED) 
#					cat("\nthe CED and the 90% confidence interval for group", 
#						gr.txt[ii], "is: \n\n", signif(sf.x * CED.matr[ii, 1], 5), 
#						"\n", signif(conf.int[ii, 1], 5), 
#						"\n", signif(conf.int[ii, 2], 5), "\n")
#			else cat("\nthe CED and the 90% confidence interval for category", 
#					CES.cat, "is: \n\n", signif(sf.x * CED.matr[1, CES.cat], 5), 
#					"\n", signif(conf.int[1], 5), 
#					"\n", signif(conf.int[2], 5), "\n")
          
        }
        
        if (profile.ans == 2) {
          
          if (dtype == 3) {
            CED.categ <- cbind(round(t(sf.x * CED.matr), 3), scores.orig[2:length(scores.orig)])
            if(track)	cat("\nCEDs (in original dose units) at (original) severity scores:\n")
            print(CED.categ)
            
          }else if (model.ans != 30 & track) {
            if (nr.aa > 1 | nr.bb > 1) 
              for (ii in 1:nr.CED) 
                cat("\nthe CED for group", gr.txt[ii], "is: \n", 
                    signif(sf.x * CED.matr[ii, 1], 5), "\n")
            else cat("\nthe CED is: \n\n", signif(sf.x * CED.matr[1, 1], 5), "\n")
          }
        }
        
        ans.scale <- 0
        
        # remove plot
        
        if (ces.ans > 1) {
          if (dtype %in% c(4, 6, 84)) {
            ans.all$CI.plt <- TRUE
            ans.all$combi.ans <- FALSE
          }
        }
        ans.all$pi.bg.matr <- pi.bg.matr
        ans.all$shift <- shift
        if (profile.ans == 1)	ans.all$conf.int <- conf.int
        
        if (track)	print("f.mm6.cat:  END")
        
        return(ans.all)
        
      })
  
}

#' Define initial parameter values for the chosen model, binary response
#' 
#' From f.start.bin of the proast61.3 package, only: tmp.quick == FALSE
#' @param ans.all list, with all results that were obtained during the analysis
#' @param adjust boolean, if TRUE the parameter start values will be adjusted;
#' default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.start.bin <- function (ans.all, adjust = FALSE, track = FALSE){
  
  if (track) print("f.start.bin")
  
  ans.all$adjust <- adjust
  ans.all$tmp.quick <- FALSE
  
  with(ans.all, {
        
        xx.tot <- x
        eps <- 1e-06
        
        if (dtype %in% c(4, 6, 84)) {
          
          xlmp <- x
          ylmp <- y
          nlmp <- nn
          fct1.lmp <- fct1
          fct2.lmp <- fct2
        }
        
        if (dtype %in% c(2, 3)) {# Laure: added in the f.start.cat function
          
          lmp.lst <- f.lump.cat(x, y, nn, fct1, fct2, dtype, twice, track = track)
          xlmp <- matrix(lmp.lst$xlmp, ncol = 1)
          ylmp <- lmp.lst$ylmp
          nlmp <- lmp.lst$nlmp
          if (length(fct1) > 1) 
            fct1.lmp <- lmp.lst$fct1.lmp
          if (length(fct2) > 1) 
            fct2.lmp <- lmp.lst$fct2.lmp
          
        }
        
        ylmp.corr <- ylmp + 0.01 * (ylmp == 0) - 0.01 * (ylmp == 1)
        
        top <- length(ylmp)
        
        if(!adjust){
          
          switch(as.character(model.ans), 
              
              # 1: null model
              '1' = {
                
                if (max(fct2) > 1) {
                  if(track)	cat("\nATTENTION: \n\nThis model is not defined for two different b parameters!")
                  ans.all$fct2 <- rep(1, length(x))
                } else {
                  aa <- sum(ylmp)/length(ylmp)
                  regr.start <- c(rep(aa, nr.aa))
                  lb <- c(rep(eps, nr.aa))
                  ub <- c(rep(1 - eps, nr.aa))
                  scale.dum <- rep(1, nr.aa)
                }
                
              }, 
              
              # 3: two-stage without CED
              '3' = {
                
                y.tmp <- logb(1 - ylmp.corr)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 1e-04
                bb <- -1/res.lm[[1]][2]
                bb <- abs(bb)
                cc <- 1
                regr.start <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                    cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), eps)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 
                    Inf)
                scale.dum <- c(rep(1, nr.aa), rep(abs(1/bb), 
                        nr.bb), abs(1/cc))
                
              },
              
              # 13: E3-CED
              '13' = {
                
                y.tmp <- logb(ylmp.corr/(1 - ylmp.corr))
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- -res.lm[[1]][1]
                bb <- res.lm[[1]][2]
                regr.start <- c(rep(aa, nr.aa), rep(bb, nr.bb))
                lb <- c(rep(-Inf, nr.aa), rep(-Inf, nr.bb))
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb))
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
                
              },
              
              # 14: full model
              '14' = {
                
                if (dtype == 2) {
                  regr.start <- as.numeric(y)
                  if(track)	cat("\nNOTE:  Full model may not be applicable for dtype = 2\n\n")
                }
                
                if (nr.aa == 1 && nr.bb == 1) covar.tmp <- F else covar.tmp <- T
                
                if (max(covariate) > 1) {
                  covar.tmp <- T
                  twice <- T
                  fct1 <- covariate
                  fct2 <- covariate
                  ans.all$fct1 <- covariate
                  ans.all$fct2 <- covariate
                }
                
                if (dtype == 6 || dtype == 4) {
                  if (!covar.tmp) {
                    x.full <- tapply(x, x, mean)
                    kk.tot <- as.numeric(tapply(kk, x, sum))
                    nn.tot <- as.numeric(tapply(nn, x, sum))
                    xx.tot <- as.numeric(tapply(x, x, sum))
                    pi.full <- kk.tot/nn.tot
                    fct1.full <- rep(1, length(x.full))
                    fct2.full <- rep(1, length(x.full))
                  }
                  if (covar.tmp) {
                    pi.full <- numeric()
                    x.full <- numeric()
                    xx.tot <- x
                    fct1.full <- numeric()
                    fct2.full <- numeric()
                    if (twice) for (jj in levels(factor(fct2))) {
                        x.tmp <- x[fct2 == jj]
                        y.tmp <- y[fct2 == jj]
                        xx.tmp <- as.numeric(tapply(x.tmp, x.tmp, mean))
                        x.full <- c(x.full, xx.tmp)
                        pi.tmp <- as.numeric(tapply(y.tmp, x.tmp, mean))
                        pi.full <- c(pi.full, pi.tmp)
                        fct1.tmp <- fct1[fct2 == jj]
                        fct1.tmp <- as.numeric(tapply(fct1.tmp, x.tmp, mean))
                        fct1.full <- c(fct1.full, fct1.tmp)
                        fct2.tmp <- fct2[fct2 == jj]
                        fct2.tmp <- as.numeric(tapply(fct2.tmp, x.tmp, mean))
                        fct2.full <- c(fct2.full, fct2.tmp)
                      }
                    if (!twice) 
                      for (jj in levels(factor(fct2))) 
                        for (ii in levels(factor(fct1))) {
                          x.tmp <- x[fct1 == ii & fct2 == jj]
                          y.tmp <- y[fct1 == ii & fct2 == jj]
                          xx.tmp <- as.numeric(tapply(x.tmp, x.tmp, mean))
                          x.full <- c(x.full, xx.tmp)
                          pi.tmp <- as.numeric(tapply(y.tmp, x.tmp, mean))
                          pi.full <- c(pi.full, pi.tmp)
                          fct1.tmp <- fct1[fct1 == ii & fct2 == jj]
                          fct1.tmp <- as.numeric(tapply(fct1.tmp, x.tmp, mean))
                          fct1.full <- c(fct1.full, fct1.tmp)
                          fct2.tmp <- fct2[fct1 == ii & fct2 == jj]
                          fct2.tmp <- as.numeric(tapply(fct2.tmp, x.tmp, mean))
                          fct2.full <- c(fct2.full, fct2.tmp)
                        }
                  }
                  regr.start <- pi.full
                }
                
                regr.start <- regr.start + eps * (regr.start == 0) - eps * (regr.start == 1)
                regr.start <- regr.start + eps * (regr.start == 0) - eps * (regr.start == 1)
                lb <- rep(eps, length(regr.start))
                ub <- rep(1 - eps, length(regr.start))
                scale.dum <- c(rep(1, length(regr.start)))
                
              }, 
              
              # E5-CED
              '15' = {
                
                y.tmp <- logb(1 - ylmp.corr)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 1e-04
                bb <- -1/res.lm[[1]][2]
                if (ces.ans == 1) CES <- 0.5
                BMD <- -bb * logb(1 - CES)
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb))
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb))
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
                
              },
              
              # 16: two-stage model
              '16' = {
                
                y.tmp <- logb(1 - ylmp.corr)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 1e-04
                bb <- -1/res.lm[[1]][2]
                bb <- abs(bb)
                cc <- 1
                
                BMD <- f.inv.bin(model.ans = 3, regr.par = c(aa, bb, cc), 
                    CES = CES, ces.ans = ces.ans, track = track)#model.ans = 3,
                
                if (ces.ans == 1) CES <- 0.5
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), eps)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 1e+12)
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), 1)
                
              }, 
              
              # 18: log-logistic model
              '18' = {
                
                y.tmp <- logb(ylmp.corr/(1 - ylmp.corr))
                dose.tmp <- xlmp
                min.dose <- dose.tmp[dose.tmp != 0]
                dose.tmp[dose.tmp == 0] <- min.dose[2]/10
                dose.tmp <- logb(dose.tmp)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 1e-04
                cc <- -abs((res.lm[[1]][2]))
                bb <- exp(res.lm[[1]][1]/cc)
                cc <- -cc
                CES.tmp <- CES/(1 - aa)
                if (ces.ans == 1) CES.tmp <- 0.5
                BMD <- bb/exp(logb((1 - CES.tmp)/CES.tmp)/cc)
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
                
              }, 
              
              # 19: Weibull model
              '19' = {
                
                y.tmp <- logb(1 - ylmp.corr)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 1e-04
                bb <- -1/res.lm[[1]][2]
                cc <- 1
                if (ces.ans == 1) CES <- 0.5
                BMD <- bb * (-logb(1 - CES))^(1/cc)
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
                
              }, 
              
              # 21: log-probit model
              '21' = {
                
                y.tmp <- qnorm(ylmp.corr)
                dose.tmp <- xlmp
                min.dose <- dose.tmp[dose.tmp != 0]
                dose.tmp[dose.tmp == 0] <- min.dose[2]/10
                dose.tmp <- logb(dose.tmp)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- ylmp[1]/nlmp[1] + 0.01
                cc <- abs(res.lm[[1]][2])
                bb <- exp(-res.lm[[1]][1]/cc)
                if (ces.ans == 1) CES <- 0.5
                BMD <- bb * exp(qnorm(CES/(1 - aa))/cc)
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
                
              }, 
              
              # H3-CED
              '23' = {
                
                y.tmp <- qnorm(ylmp.corr)
                dose.tmp <- xlmp
                min.dose <- dose.tmp[dose.tmp != 0]
                dose.tmp[dose.tmp == 0] <- min.dose[2]/10
                dose.tmp <- logb(dose.tmp)
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                bb <- res.lm[[1]][2]
                aa <- exp(-res.lm[[1]][1]/bb)
                BMD <- logb(aa) + qnorm(CES)/bb
                BMD <- exp(BMD)
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb))
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb))
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
                
              },
              
              # 24: Gamma model
              '24' = {
                
                BMD <- mean(x)
                aa <- 0.01
                cc <- 1
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
                lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
                
              }, 
              
              # 25: probit model for quantal, H5-CED for ordinal
              '25' = {
                
                y.tmp <- qnorm(ylmp.corr)
                dose.tmp <- xlmp
                min.dose <- dose.tmp[dose.tmp != 0]
                dose.tmp[dose.tmp == 0] <- min.dose[2]/10
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                # Proast 62.10
                aa <- res.lm[[1]][1]
                bb <- res.lm[[1]][2]
                # Proast 61.3
#                bb <- res.lm[[1]][2]
#                aa <- (-res.lm[[1]][1]/bb)
                if (ces.ans == 1) {
                  BMD <- mean(xlmp)
                  # Proast 62.10
                  BMD <- abs(-aa/bb)
#                  regr.start <- c(rep(BMD, nr.aa), rep(bb, nr.bb))
                  # Proast 62.10
                  regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
                  lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
                  ub <- c(rep(Inf, nr.aa), rep(Inf, nr.bb))
                  scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
                }
                if (ces.ans %in% 2:3) {
                  BMD <- qnorm(CES * (1 - pnorm(-aa * bb)) + pnorm(-aa * bb))
#                  BMD <- BMD/bb + aa
                  # Proast 62.10
                  BMD <- abs(BMD/bb + aa)
                  BMD <- BMD * 2
                  nr.aa <- 1
                  nr.bb <- 1
                  regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
                  lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
                  ub <- c(rep(Inf, nr.aa), rep(Inf, nr.bb))
                  scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
                }
                
              }, 
              
              # 26: logistic model
              '26' = {
                
                y.tmp <- logb(ylmp.corr/(1 - ylmp.corr))
                dose.tmp <- xlmp
                dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
                res.lm <- lm(yyy ~ dose.tmp, dt.fr)
                aa <- -res.lm[[1]][2]
                bb <- res.lm[[1]][2]
                if (ces.ans == 1) CES <- 0.5
                BMD <- -logb((1 - CES)/(CES + exp(aa))) - aa
                BMD <- BMD/bb
                regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
                lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
                ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb))
                scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
              }
          
          )
          
          par.start <- regr.start
          
          if (dtype == 6) {
            if (!(model.ans == 14 && model.type == 1))
              par.start <- c(alfa.start, regr.start)
            else par.start <- c(10, regr.start)
            lb <- c(1e-10, lb)
            ub <- c(Inf, ub)
          }
          
          nrp <- length(regr.start)
          
          
          
          # cens parameter removed
          loglik.old <- -f.lik.cat(theta = par.start, 
              x = x, y = y, kk = kk, 
              nn = nn, dtype = dtype, 
              fct1 = fct1, fct2 = fct2, 
              nrp = nrp, nth = 1, 
              nr.aa = nr.aa, nr.bb, 
              # Note: in original code, model.type <- 1
              model.ans = model.ans,  model.type = model.type, CES = CES, 
              ttt = ttt, twice = twice, alfa.length = alfa.length, 
              x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
              ces.ans = ces.ans, CES1 = CES1, CES2 = CES2, 
              nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot,
              track = track)
          
          if(track)	cat("\nLog-likelihood value at start values: ", signif(loglik.old, 5), "\n")
          if(track)	cat("\nModel has not yet been fitted !\n")
          ans.all$loglik.old <- loglik.old
          
        } else { # if adjust == TRUE
          
          loglik.old <- 0
          ans <- 1
          ans.all.tmp <- ans.all
          ans.all.tmp$x <- xlmp
          ans.all.tmp$y <- ylmp
          ans.all.tmp$nn <- nlmp
          ans.all.tmp$heading <- "starting values"
          ans.all.tmp$l.ty <- 2
          ans.all.tmp$xy.lim[3] <- max(xlmp)
          
          par.start <- startValues
          ans.all.tmp$regr.par <- par.start
          if (dtype == 6) 
            ans.all.tmp$regr.par <- par.start[-1]
          
          loglik.new <- NA
          
          loglik.new <- -f.lik.cat(par.start, x, y, kk, 
              nn, dtype, fct1, fct2, nrp, nth = 1, nr.aa, 
              nr.bb, model.ans, model.type = 1, CES = CES, 
              ttt = ttt, twice = twice, alfa.length = alfa.length, 
              x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
              ces.ans = ces.ans, CES1 = CES1, CES2 = CES2, 
              nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot)
          
          if(is.na(loglik.new))
            ans.all$errorAdjustStartValues <- TRUE
#            stop("Please choose other start values for the model parameters:\n The log-likelihood could not be calculated")
          
          ans.all.tmp$show <- ""
          ans.all.tmp$CED.matr <- NA
          loglik.old <- loglik.new
          ans.all$fitted <- FALSE          
          
        }
        
        ans.all$par.start <- par.start
        ans.all$regr.start <- regr.start
        ans.all$nrp <- nrp
        ans.all$npar <- length(ans.all$par.start)
        
        if (model.ans == 14 & dtype %in% c(4, 6)) {
          
          ans.all$x.full <- x.full
          ans.all$fct1.full <- fct1.full
          ans.all$fct2.full <- fct2.full
          
        }
        
        if(!adjust){
          
          ans.all$lb <- lb
          ans.all$ub <- ub
          if (dtype == 6) 
            betabin <- 2
          else betabin <- 1
          if (model.ans %in% c(18:19, 24)) 
            if (constr != -Inf) 
              ans.all$lb[nr.aa + nr.bb + betabin] <- constr
          scale.dum <- abs(1/par.start)
          scale.dum[scale.dum < 0.001] <- 0.001
          scale.dum[scale.dum > 1000] <- 1000
          ans.all$scale.dum <- scale.dum
          ans.all$nn.tot <- nn.tot
          ans.all$kk.tot <- kk.tot
          ans.all$xx.tot <- xx.tot
          
        }
        
        
        if (track)	print("f.start.bin:  END")
        
        return(ans.all)
        
      })
}

#' ? Add values to e.g. response y to avoid log of 0
#' @param x vector independent variable
#' @param y vector response
#' @param nn vector number of observations
#' @param fct1 vector covariate on which parameter a is dependent
#' @param fct2 vector covariate on which parameter b is dependent
#' @param dtype integer response data type
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with some values
#' @export
f.lump.cat <- function (x, y, nn, fct1, fct2, dtype, twice, track = FALSE) {
  
  if (track)	print("f.lump.cat")
  
  x.ss <- tapply(x, x, length)
  if (mean(x.ss) > 8) 
    ngroups <- length(x.ss)
  else ngroups <- 5
  xlmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
  ylmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
  nlmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
  fct1.lmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
  fct2.lmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
  kk <- 0
  for (jj in 1:max(fct2)) {
    if (!twice) {
      for (ii in 1:max(fct1)) {
        kk <- kk + 1
        x.part <- x[fct1 == ii & fct2 == jj]
        y.part <- y[fct1 == ii & fct2 == jj]
        sub <- round(length(x.part)/ngroups)
        e1 <- 0
        e2 <- 0
        for (ee in 1:ngroups) {
          e2 <- e2 + sub
          xtmp <- x.part[e1:e2]
          ytmp <- y.part[e1:e2]
          xtmp <- xtmp[is.finite(xtmp)]
          ytmp <- ytmp[is.finite(ytmp)]
          xlmp[kk, ee] <- mean(xtmp)
          ylmp[kk, ee] <- mean(ytmp)
          e1 <- e1 + sub
        }
        nlmp[kk, ] <- rep(sub, ngroups)
        fct1.lmp[kk, ] <- rep(ii, ngroups)
        fct2.lmp[kk, ] <- rep(jj, ngroups)
      }
    }
    else {
      kk <- kk + 1
      x.part <- x[fct2 == jj]
      y.part <- y[fct2 == jj]
      sub <- round(length(x.part)/ngroups)
      e1 <- 0
      e2 <- 0
      for (ee in 1:ngroups) {
        e2 <- e2 + sub
        xtmp <- x.part[e1:e2]
        ytmp <- y.part[e1:e2]
        xtmp <- xtmp[is.finite(xtmp)]
        ytmp <- ytmp[is.finite(ytmp)]
        xlmp[kk, ee] <- mean(xtmp)
        ylmp[kk, ee] <- mean(ytmp)
        e1 <- e1 + sub
      }
      fct1.lmp[kk, ] <- rep(jj, ngroups)
      fct2.lmp[kk, ] <- rep(jj, ngroups)
      nlmp[kk, ] <- rep(sub, ngroups)
    }
  }
  xlmp <- matrix(t(xlmp), ncol = 1, nrow = ngroups * kk)
  ylmp <- matrix(t(ylmp), ncol = 1, nrow = ngroups * kk)
  nlmp <- matrix(t(nlmp), ncol = 1, nrow = ngroups * kk)
  fct1.lmp <- matrix(t(fct1.lmp), ncol = 1, nrow = ngroups * kk)
  fct2.lmp <- matrix(t(fct2.lmp), ncol = 1, nrow = ngroups * kk)
  ylmp <- ylmp[!is.na(xlmp)]
  nlmp <- nlmp[!is.na(xlmp)]
  fct1.lmp <- fct1.lmp[!is.na(xlmp)]
  fct2.lmp <- fct2.lmp[!is.na(xlmp)]
  xlmp <- xlmp[!is.na(xlmp)]
  lmp.lst <- list(xlmp = xlmp, ylmp = ylmp, nlmp = nlmp, fct1.lmp = fct1.lmp, 
      fct2.lmp = fct2.lmp, dtype = dtype)
  
  if (track)	print("f.lump.cat end: END")
  
  return(lmp.lst)
  
}

#' return the CED value for a set of parameters, CES and ces.ans
#' 
#' Originally the function can deal all model.ans, 
#' but because currently this function is only used:
#' \itemize{
#' \item{in the \code{f.start.bin} function where model.ans == 16, in which model.ans is set to 3}
#' \item{in \code{f.ced.cat} for the categorical models: 
#' 	\itemize{
#' 	\item{quantal/binary response: }{(1, 14), 16, 18, 19, 21, 24, 25, 26}
#'  \item{ordinal response: }{13, 15, 23, 25 so only ces.ans == 1}
#'  }}
#' } 
#' Only this part of the code is retained. 
#' @param model.ans integer, type of response model, only c(3 available
#' @param regr.par vector with parameters values (a, b, c)
#' @param CES numeric, CES value
#' @param ces.ans integer, type of benchmark response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return CED value
#' @export
f.inv.bin <- function (model.ans, regr.par, CES, ces.ans, track = FALSE){
  
  if (track)	print("f.inv.bin")
  
  aa <- regr.par[1]
  bb <- regr.par[2]
  cc <- regr.par[3]
  dd <- regr.par[4]
  
  if(model.ans == 1){
    if(track)	cat("BMD calculation not implemented for null model\n")
    CED <- NA
  }
  
  if (model.ans == 14) {
    if(track)	cat("BMD calculation not implemented for full model\n")
    CED <- NA
  }
  
  if (model.ans %in% 15:28)	CED <- bb
  
  # Proast version 61.3
  # MV changed this into bb for ces.ans > 1
#  if (model.ans == 25)
#    if(ces.ans == 1) 
#      CED <- aa else
#      CED <- bb
  # Proast version 62.10
  if (model.ans == 25) 
    CED <- bb
  
  if (ces.ans == 1){
    
    switch(as.character(model.ans),
        
        '3' = {
          
          qq <- log(2 - 2 * aa)
          CED <- sqrt((bb^2/cc) * (qq + 1/(4 * cc)))
          CED <- CED - 0.5 * bb/cc
          
        },
        
        '13' = CED <- -aa/bb
    
    )
    
  }
  
  
  if(model.ans == 3){
    
    if (ces.ans == 2){
      
      # only code for model.ans == 3
      if (bb == 0){
        if(track)	cat("\n\nATTENTION: (one of the) parameter(s) bb  is zero!\n\n \nTherefore the CED cannot be calculated\n") 
      }else if(cc == 0){
        CED <- -bb * logb(1 - CES/(1 - aa)) 
      }else {
        dum1 <- (-0.5 * bb)/cc
        dum2 <- sqrt((0.25 * bb^2)/cc^2 - (bb^2/cc) * logb(1 - CES/(1 - aa)))
        CEDa <- dum1 + dum2
        CEDb <- dum1 - dum2
        CED <- max(CEDa, CEDb)
      }
    }
    
    if (ces.ans == 3){
      
      # only code for model.ans == 3	
      if (bb == 0){
        if(track)	cat("\n\nATTENTION: (one of the) parameter(s) bb is zero!\n\nTherefore the CED cannot be calculated\n") 
      }else if (cc < 1e-05){
        CED <- -bb * logb(1 - CES) 
      }else {
        dum1 <- (-0.5 * bb)/cc
        dum2 <- sqrt((0.25 * bb^2)/cc^2 - (bb^2/cc) * logb(1 - CES))
        CEDa <- dum1 + dum2
        CEDb <- dum1 - dum2
        CED <- max(CEDa, CEDb)
      }
    }	
    
  }
  
  if (track)	print("f.inv.bin: END")
  
  return(CED)
  
}


#' Calculate likelihood for categorical response
#' @param theta numeric vector, the initial regression parameter values
#' @param x numeric vector, the dose values
#' @param y numeric vector, the response values 
#' @param kk vector of parameters?
#' @param nn vector of total number
#' @param dtype integer, determines the type of response
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param nrp values?
#' @param nth values?
#' @param nr.aa values?
#' @param nr.bb values?
#' @param model.ans integer, type of response model, only 3 available
#' @param model.type type of model
#' @param CES numeric, value for the CES
#' @param ttt numeric, time variable 
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param alfa.length length of alpha
#' @param x.full numeric vector, the dose values
#' @param fct1.full numeric, value for parameter a
#' @param fct2.full numeric, value for parameter b
#' @param ces.ans index type of benchmark response
#' @param fct3 numeric, value for parameter c
#' @param fct4 numeric, value for parameter d
#' @param fct5 numeric, value for parameter e
#' @param CES.cat value
#' @param decr.zz is z decreasing?
#' @param CES1 numeric, value for the first CES
#' @param CES2 numeric, value for the second CES
#' @param nn.tot vector of total number
#' @param kk.tot vector of parameters?
#' @param xx.tot vector of concentrations
#' @param fct3.ref reference for parameter 3
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric value, minus the sum of the scores (total likelihood)
#' @export
f.lik.cat <- function (theta, x, y = NA, kk, nn, dtype, fct1 = 1, fct2 = 1, 
    nrp, nth = NA, nr.aa = 1, nr.bb = 1, model.ans, model.type,
    CES = NA, ttt = 0, twice = NA, 
    alfa.length = 0, x.full = NA, 
    fct1.full = NA, fct2.full = NA, ces.ans = 1, 
    fct3 = 1, fct4 = 1, fct5 = 1, CES.cat = 1,
    decr.zz = TRUE, 
    CES1 = CES1, CES2 = CES2, 
    nn.tot, kk.tot, xx.tot = xx.tot, fct3.ref, track = FALSE){
  
  if (track)	print("f.lik.cat")
  
  if (sum(is.na(theta)) > 0)	return(NA)
  
  if (model.type == 1 & model.ans != 25 & sum(theta == 0, na.rm = T) > 0) {
    theta[theta == 0] <- theta[theta == 0] + 0.1
  }
  
  if (sum(is.na(theta) > 0) > 0) {
    theta[is.na(theta)] <- 0.1
  }
  
  par.lst <- f.split.par(MLE = theta, nrp, nth, dtype, track = track)
  regr.par <- par.lst$regr.par
  th.par <- par.lst$th.par
  sig.par <- par.lst$sig.par
  
  
  if (dtype != 6) {
    
    # Added one line in new proast
    if(!(model.ans %in% c(12,13,25,26)))   # for probit and logit model parameter a can be zero
      if (model.type == 1 & prod(theta[1:nr.aa]) == 0) {
        kk <- kk[x != 0]
        nn <- nn[x != 0]
        fct1 <- fct1[x != 0]
        fct2 <- fct2[x != 0]
        x <- x[x != 0]
      }
    
    if (model.type == 1 & model.ans == 14 && nr.aa == 1 && nr.bb == 1) {
      kk <- kk.tot
      nn <- nn.tot
      x <- xx.tot
    }
    
    if(model.type == 1){
      
      pipi <- f.expect.bin(model.ans = model.ans, x = x, 
          regr.par = theta, fct1 = fct1, fct2 = fct2, 
          #			name = F, 
          kk = kk, nn = nn, dtype = dtype, 
          CES = CES, ttt = ttt, twice = twice, 
          x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
          #			trace = FALSE, 
          ces.ans = ces.ans, CES1 = CES1, 
          CES2 = CES2, track = track)
      
      if (sum(!(is.na(pipi) == 0)))	return(1e+10)
      
    }
    
    if (model.type == 2) {
      
      pipi <- f.expect.cat(model.type, model.ans, x, regr.par, 
          th.par, sig.par, fct1, fct2, ttt = ttt, twice = twice, 
          fct3 = fct3, ces.ans = ces.ans, CES = CES, CES.cat = CES.cat, 
          decr.zz = decr.zz, dtype = dtype, fct4 = fct4, 
          fct5 = fct5, fct3.ref = fct3.ref,
          track = track)
      
    }
    
    
    lik <- numeric(0)
    
    if (dtype == 4 | dtype == 6 | dtype == 84) {
      
      if (sum(is.na(pipi)) > 0) {
        
        totlik <- -1e+12
        
      }else {
        
        if ((min(pipi) == 0) | (max(pipi) == 1)) {
          nth <- 1
          pi.tmp <- rep(pipi[1], nn[1])
          y <- c(rep(1, kk[1]), rep(0, nn[1] - kk[1]))
          for (ii in 2:length(x)) {
            pi.tmp <- c(pi.tmp, rep(pipi[ii], nn[ii]))
            y <- c(y, rep(1, kk[ii]), rep(0, nn[ii] - kk[ii]))
          }
          pipi <- matrix(pi.tmp, ncol = 1)
          dtype <- 400
          
        }else {
          
          if (max(pipi) > 1) {
            if(track)	cat("\n\nATTENTION: pi > 1 \n")
            if(track)	cat("current parameter values:\n")
            print(signif(theta, 3))
          }
          
          # sum(cens) == 0)
          loglik <- kk * logb(pipi) + (nn - kk) * logb(1 - pipi)
          totlik <- sum(loglik)
          
#					if (sum(cens) > 0) {
#						lik <- dbinom(kk, nn, pipi) * (cens != 1) + 
#							(1 - pbinom(kk, nn, pipi) + 
#							dbinom(kk, nn, pipi)) * (cens == 1)
#						totlik <- sum(logb(lik))
#					}
          
        }
        
      }
      
    }
    
    if (dtype == 2) {# Laure: added for dtype == 3
      nth <- 1
      pipi <- matrix(pipi, ncol = 1)
    }
    
    if (dtype == 2 | dtype == 3 | dtype == 400) {
      
      pi0 <- rep(1, length(y))
      for (k in 1:nth) pi0 <- pi0 - pipi[, k]
      pi0[pi0 == 0] <- 1e-12
      lik <- (y == 0) * pi0
      for (k in 1:nth) lik <- lik + (y == k) * pipi[, k]
      loglik <- logb(lik)
      totlik <- sum(loglik)
      
    }
    
    if (!is.na(totlik)) {
      if (totlik == -Inf) {
        totlik <- -2e+10
      }
    }
  }
  
  if (dtype == 6) {
    
    alfa <- theta[1:alfa.length]
    
    if (alfa.length > 1) {
      xtot <- as.numeric(tapply(x, x, mean))
      xtot <- round(xtot, 10)
      x.tmp <- round(x, 10)
      alfa.tmp <- rep(0, length(x.tmp))
      for (ii in 1:length(xtot)) 
        alfa.tmp <- alfa.tmp + alfa[ii] * (x.tmp == xtot[ii])
      alfa <- alfa.tmp
    }
    
    if (model.type == 1)
      pipi <- f.expect.bin(model.ans = model.ans, x = x, regr.par = regr.par, 
          fct1 = fct1, fct2 = fct2, 
          #			name = F, 
          kk = kk, nn = nn, dtype = dtype, 
          CES = CES, ttt = ttt, 
          twice = twice, x.full = x.full, 
          fct1.full = fct1.full, fct2.full = fct2.full, 
          #trace = trace, 
          ces.ans = ces.ans, CES1 = CES1, CES2 = CES2,
          track = track)
    
    if (model.type == 2) 
      pipi <- f.expect.cat(model.type, model.ans, x, regr.par, 
          th.par, sig.par, fct1, fct2, twice = twice, CES = CES, 
          CES.cat = CES.cat, fct3 = 1, ces.ans = ces.ans, 
          dtype = dtype, fct4 = fct4, fct5 = fct5, fct3.ref = fct3.ref, 
          track = track)
    
    pipi[pipi < 1e-04] <- 0.001
    pipi[pipi > 0.9999] <- 1 - 0.001
    lik <- numeric(0)
    beta.0 <- (alfa * (1 - pipi))/pipi
    log.beta.numerator <- lbeta(kk + alfa, nn - kk + beta.0)
    log.beta.denominator <- lbeta(alfa, beta.0)
    bin.coef <- 1
    loglik <- log(bin.coef) + log.beta.numerator - log.beta.denominator
    totlik <- sum(loglik)
    if (!is.finite(totlik)) 
      totlik <- NA

  }
  
  if (track)	print("f.lik.cat: END")
  
  return(-totlik)
  
}


#' split the vector of parameters
#' @param MLE vector of parameters
#' @param nrp values?
#' @param nth values?
#' @param dtype integer, determines the type of response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with parameters split by category
#' @export
f.split.par <- function (MLE, nrp, nth, dtype, track = FALSE){
  
  if (track)	print("f.split.par")
  
  if (dtype == 6) {
    MLE <- MLE[-1]
  }
  
  regr.par <- MLE[1:nrp]
  th.par <- MLE[(nrp + 1):(nrp + nth)]
  sig.par <- MLE[nrp + nth + 1]
  
  if (track)	print("f.split.par: END")
  
  return(list(regr.par = regr.par, th.par = th.par, sig.par = sig.par))
  
}

#' Calculate expected response values, for categorical response
#' @param model.ans integer, determines the type of model to be fitted
#' @param x numeric vector, the dose values
#' @param regr.par numeric vector, regression parameter values
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param kk vector of parameters?
#' @param nn vector of total number
#' @param dtype integer, determines the type of response
#' @param CES numeric, value for the CES
#' @param ttt numeric, time variable (redundant)
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param x.full vector of concentrations
#' @param fct1.full numeric, value for parameter a
#' @param fct2.full numeric, value for parameter b
#' @param ces.ans index type of benchmark response
#' @param CES1 numeric, value for the first CES (redundant)
#' @param CES2 numeric, value for the second CES (redundant)
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, the expected response values under the estimated 
#' regression model
#' @export
f.expect.bin <- function (model.ans = NA, x, regr.par = 0, fct1 = 1, fct2 = 1, 
    kk = -1, nn = -1, dtype = 6, CES = NA, ttt = 0, 
    twice = F, x.full = NA, fct1.full = NA, fct2.full = NA, #trace = F, 
    ces.ans = 3, CES1, CES2, track = FALSE) {
  
  if (track)	print("f.expect.bin")
  
  nr.aa <- max(fct1)
  nr.bb <- max(fct2)
  nrp <- length(regr.par)
  aa0 <- rep(0, length(x))
  bb0 <- rep(0, length(x))
  aa.tmp <- regr.par[1:nr.aa]
  bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
  
  for (ii in 1:nr.aa) 
    aa0 <- aa0 + aa.tmp[ii] * (fct1 == ii)
  
  for (jj in 1:nr.bb)
    bb0 <- bb0 + bb.tmp[jj] * (fct2 == jj)
  
  cc <- regr.par[nr.aa + nr.bb + 1]
  dd <- regr.par[nr.aa + nr.bb + 2]
  ee <- regr.par[nr.aa + nr.bb + 3]
  
  switch(as.character(model.ans), 
      
      # 1: null model
      '1' = pipi <- aa0, 
      
      # 3: two-stage model (without bmd)
      '3' = pipi <- aa0 + (1 - aa0) * (1 - exp(-x/bb0 - cc * (x/bb0)^2)),
      
      # 13: E3-CED
      '13' = pipi <- 1/(1 + exp(-aa0 - bb0 * x)),
      
      # 14: full model
      '14' = {
        
        if (dtype == 6) {
          
          pipi <- rep(0, length(x))
          
          if (twice) {
            for (ii in (1:max(fct1.full))) {
              regr.par.tmp <- regr.par[fct1.full == ii]
              x.full.tmp <- x.full[fct1.full == ii]
              for (kk in 1:length(x.full.tmp)) 
                pipi <- pipi + regr.par.tmp[kk] * (x == x.full.tmp[kk]) * (fct1 == ii)
            }
          }else{
            for (jj in (1:max(fct2.full))) 
              for (ii in (1:max(fct1.full))) {
                regr.par.tmp <- regr.par[fct1.full == ii & fct2.full == jj]
                x.full.tmp <- x.full[fct1.full == ii & fct2.full == jj]
                for (kk in 1:length(x.full.tmp))
                  pipi <- pipi + regr.par.tmp[kk] * 
                      (x == x.full.tmp[kk]) * (fct1 == ii) * (fct2 == jj)
              }
          }
        }
        
        if (dtype %in% c(2, 4))	pipi <- regr.par # LAURE: dtype == 2 added because no code in the original function 
        
      }, 
      
      # 15: E5-CED
      '15' = {
        
        if (is.na(CES)) cat("\nValue for CES is not known !! \n")
        if (ces.ans == -1) 
          pipi <- aa0 + (1 - aa0) * (1 - exp((x * logb(1 - CES))/bb0)) else
        {
          dum0 <- f.bb.bin(model.ans, aa0, cc = NA, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
          pipi <- aa0 + (1 - aa0) * (1 - exp(-x/dum0))
        }
        
      },
      
      # 16: two-stage model
      '16' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        
        if (ces.ans == -1) {
          if (cc == 0) {
            pipi <- aa0 + (1 - aa0) * (1 - exp((x * logb(1 - CES))/bb0))
          } else {
            dum0 <- ((-2 * bb0 * cc)/(1 - sqrt(1 - 4 * cc * logb(1 - CES))))
            pipi <- aa0 + (1 - aa0) * (1 - exp(-x/dum0 - cc * (x/dum0)^2))
          }
        } else {
          dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
          pipi <- aa0 + (1 - aa0) * (1 - exp(-x/dum0 - cc * (x/dum0)^2))          
        }
        
      }, 
      
      # 18: logistic model
      '18' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        if (ces.ans == -1)
          dum0 <- bb0 * exp(logb((1 - CES)/CES)/cc) 
        else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
        
        pipi <- aa0 + (1 - aa0) * (1/(1 + exp(cc * (logb(dum0) - logb(x)))))
        
      }, 
      
      # 19: Weibull model
      '19' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        if (ces.ans == -1)
          dum0 <- bb0/(-logb(1 - CES))^(1/cc)
        else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
        
        pipi <- aa0 + (1 - aa0) * (1 - exp(-(x/dum0)^cc))
        
      }, 
      
      # 21: log-probit model
      '21' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        if (ces.ans == -1)
          dum0 <- bb0 * exp(-qnorm(CES)/cc)
        else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
        
        pipi <- aa0 + (1 - aa0) * pnorm((logb(x) - logb(dum0)) * cc)
        
      }, 
      
      # 23: H3-CED
      '23' = {
        
        if (is.na(CES)) cat("\nValue for CES is not known !! \n")
        dum0 <- qnorm(CES)/(logb(bb0) - logb(aa0))
        pipi <- pnorm((logb(x) - logb(aa0)) * dum0)
        
      },
      
      # 24: Gammma model
      '24' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        if (ces.ans == -1)
          dum0 <- qgamma(CES, cc)/bb0
        else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
              CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
        
        pipi <- aa0 + (1 - aa0) * pgamma(dum0 * x, cc)
        
      }, 
      
      # 25: probit model for quantal/binary, H5-CED for ordinal
      '25' = {
        
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        
        if (ces.ans > 1 & ((nr.aa > 1) | (nr.bb > 1))) {
          
          pipi <- NA
          stop("Factor for parameter 'a' or 'b' not implemented for probit model")
          ans.all$fct1 <- 1
          ans.all$fct2 <- 1
          
        } else {
          
          aa <- aa0[1]
          bb <- bb0[1]
          # Proast version 62.10
          dum <- f.bb.bin(model.ans, aa, cc = NA, dd = NA, 
              CED = bb, CES = CES, ces.ans = ces.ans, track = track)
          pipi <- pnorm(aa + dum * x)
          # Proast version 61.3
#          if (ces.ans %in% 2:3) {
#            dum <- f.bb.bin(model.ans, aa, cc = NA, dd = NA, 
#                CED = bb, CES = CES, ces.ans = ces.ans, track = track)
#            pipi <- pnorm(((x) - (aa)) * dum)
#          }
#          
#          if (ces.ans == 1) {
#            dum0 <- f.bb.bin(model.ans, aa, cc = cc, 
#                dd = NA, CED = bb, CES = CES, ces.ans = ces.ans, track = track)
#            pipi <- pnorm(((x) - (dum0)) * bb0)
#          }
          
        }
        
      }, 
      
      # 26: logistic model
      '26' = {
        if (is.na(CES) & track) cat("\nValue for CES is not known !! \n")
        dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
            CED = bb0, CES = CES, ces.ans = ces.ans, track = track)
        pipi <- 1/(1 + exp(-aa0 - dum0 * x))
      }
  
  )
  
  if (track)	print("f.expect.bin: END")
  
  
  return(pipi)
  
}


#' compute the value for parameter b, from values of other parameters
#' 
#' From the original function from proast61.3, only the code for 
#' model.ans = [15, 16, 18, 19, 21, 24, 25, 26] is retained.
#' @param model.ans integer, determines the type of model to be fitted
#' @param aa value for parameter a
#' @param cc value for parameter c
#' @param dd value for parameter d
#' @param CED numeric, value for the CED
#' @param CES numeric, value for the CES
#' @param ces.ans index type of benchmark response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return value for the b parameter 
#' @export
f.bb.bin <- function (model.ans, aa, cc, dd, CED, CES, ces.ans, track = FALSE) {
  
  if (track)	print("f.bb.bin")
  
  # model.ans = [15, 16, 18, 19, 21, 24, 25, 26]
  # -> [1, 2, 4, 5, 7, 10, 11, 12]
  model.ans <- model.ans - 14
  
  # Proast 62.10
  if (model.ans == 25 - 14) 
    bb <- aa
  
  if (ces.ans == 5 && model.ans == 4) {
    bb <- CED * ((1/aa - CES - 1)/CES)^(1/cc)
  }
  
  if (ces.ans == 1) 
    
    switch(as.character(model.ans),
        
        # 1
        '1' = bb <- CED/log(2 - 2 * aa),	
        
        # 2
        '2' = {
          qq <- log(2 - 2 * aa)
          bb <- sqrt((CED^2/qq) * (cc + 1/(4 * qq)))
          bb <- bb + 0.5 * CED/qq
        }, 
        
        # 4
        '4' = bb <- CED/(1 - 2 * aa)^(1/cc),
        
        # 5
        '5' = bb <- CED/(log(2 - 2 * aa))^(1/cc), 
        
        # 7	
        '7' = bb <- CED/exp(qnorm((0.5 - aa)/(1 - aa))/cc), 
        
        # 10
        '10' = {
          bb <- qgamma((0.5 - aa)/(1 - aa), cc)
          bb <- bb/CED
        }, 
        
        # 11
      # Proast version 61.3
#        '11' = bb <- aa, 
      # Proast version 62.10
        '11' = bb <- -aa/CED,
        # 12
        '12' = bb <- -aa/CED
    
    )
  
  if (ces.ans == 2) 
    
    switch(as.character(model.ans), 
        
        # 1
        '1' = {
          
          if (sum(CED == 0) > 0) {
            if(track)	cat("\n\nATTENTION: one of the CEDs is zero!\n\nTherefore bb cannot be calculated\n")
            bb <- NA
          } else {
            qq <- log(1 - CES/(1 - aa))
            bb <- -CED/qq
          }
          
        },
        
        # 2	
        '2' = {
          qq <- log(1 - CES/(1 - aa))
          if (sum(CED == 0) > 0) {
            if(track)	cat("\n\nATTENTION: one of the CEDs is zero!\n\nTherefore bb cannot be calculated\n")
            bb <- NA
          } else if (cc == 0) bb <- -CED/qq 
          else {
            # MV changed this: mistake!
            bb <- CED * sqrt(-cc/qq + 1/(4 * qq * qq))
#            bb <- CED * sqrt(cc/qq + 1/(4 * qq * qq))
            bb <- bb - CED/(2 * qq)
          }
        }, 
        
        # 4
        '4' = bb <- CED * exp(logb((1 - aa)/CES - 1)/cc), 
        
        # 5
        '5' = bb <- CED/(-logb(1 - CES/(1 - aa)))^(1/cc), 
        
        # 7
        '7' = bb <- CED/exp(qnorm(CES/(1 - aa))/cc), 
        
        # 10
        '10' = bb <- qgamma(CES/(1 - aa), cc)/CED, 
        
        # 11	
        '11' = bb <- f.uniroot.probit(aa, CED, CES, ces.ans, track = track), 
        
        # 12	
        '12' = {
          bb <- (1 + exp(-aa))/(1 + CES + CES * exp(-aa))
          bb <- bb - 1
          bb <- -logb(bb) - aa
          bb <- bb/CED
        }
    
    )
  
  
  if (ces.ans == 3) 
    
    switch(as.character(model.ans), 
        
        # 1
        '1' =  bb <- -CED/logb(1 - CES),
        
        # 2	
        '2' = {
          if (sum(CED == 0) > 0) {
            if(track)	cat("\n\nATTENTION: (one of the) CED(s) is zero!\n\nTherefore bb cannot be calculated\n")
            bb <- NA
          } else if (cc == 0) bb <- -CED/logb(1 - CES)
          else {
            bb <- ((-2 * CED * cc)/(1 - sqrt(1 - 4 * cc * logb(1 - CES))))
          }
        }, 
        
        # 4	
        '4' = bb <- CED * exp(logb((1 - CES)/CES)/cc), 
        
        # 5	
        '5' = bb <- CED/((-logb(1 - CES))^(1/cc)), 
        
        # 7	
        '7' = bb <- CED/exp(qnorm(CES)/cc), 
        
        # 10	
        '10' = bb <- qgamma(CES, cc)/CED,
        
        # 11	
        '11' = bb <- f.uniroot.probit(aa, CED, CES, ces.ans, track = track), 
        
        # 12 
        '12' = bb <- (-logb((1 - CES)/(CES + exp(aa))) - aa)/CED
    
    )
  
  if (any(is.infinite(bb)) & track) {
    cat("\n\nATTENTION: parameter values lead to infinite value for parameter b\n",
        "Results will not be reliable\n\n")
  }
  
  return(bb)
  
}

#' compute the value for the parameter b for a probit function
#'
#' From the original function from proast61.3, only the code for 
#' model.ans = 11 (original model.ans == 25), used in f.bb.bin is retained.
#' @param aa value for parameter a
#' @param CED value for the CES
#' @param CES value for the CED
#' @param ces.ans index type of benchmark response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return value for parameter b for probit function
#' @export
f.uniroot.probit <- function (aa, CED, CES, ces.ans, track = FALSE) {
  if (exists("track2")) 
    print("f.uniroot.probit")
  if (CED <= 0) {
    cat("ATTENTION (f.uniroot):  BMD smaller than (or equal to) zero\n")
    return(NA)
  }
  f.b <- function(bb, CED, aa, CES, ces.ans) {
# Proast 61.3
#    if (ces.ans == 2) 
#      zz <- (1/bb) * qnorm(CES + pnorm(-aa * bb)) + aa/bb
#    if (ces.ans == 3) 
#      zz <- (1/bb) * qnorm(CES * (1 - pnorm(-aa * bb)) + 
#              pnorm(-aa * bb)) + aa/bb
#    return(zz - CED)
    # Proast 62.10
    if (ces.ans == 2) 
      CED.appr <- (qnorm(CES + pnorm(aa)) - aa)/bb
    if (ces.ans == 3) 
      CED.appr <- (qnorm(CES * (1 - pnorm(aa)) + pnorm(aa)) - 
            aa)/bb
    return(CED.appr - CED)
  }
  # Proast version 61.3
#  bb.low <- 0.1
  # Proast 62.10
  bb.low <- 1e-05
  bb.upp <- 10
#  if (model.ans == 25) {
#    bb.low <- -10
#    bb.upp <- 10
#  }
  try1 <- f.b(bb.low, CED, aa, CES, ces.ans)
  ii <- 1
  while (is.na(try1) & ii < 10) {
    # Proast version 61.3
#    bb.low <- -(6 - ii)/aa
    # Proast 62.10
    bb.low <- bb.low/10
    try1 <- f.b(bb.low, CED, aa, CES, ces.ans)
    ii <- ii + 1
  }
  try2 <- f.b(bb.upp, CED, aa, CES, ces.ans)
  # Proast 61.3
#  if (is.na(try1)) 
#    return(NA)
#  if (is.na(try2)) 
#    return(NA)
# Proast 62.10
while (is.na(try2) & ii < 10) {
  bb.upp <- bb.upp * 10
  try1 <- f.b(bb.upp, CED, aa, CES, ces.ans)
  ii <- ii + 1
}
  if (sign(try1) == sign(try2)) {
    # Proast 61.3
#    cat("   (f.uniroot.probit):  no root found for parameter b\n")
  # Proast 62.10
  cat("   (f.uniroot.probit):  no root found for parameter a\n")
    return(NA)
  }
  root.out <- uniroot(f.b, interval = c(bb.low, bb.upp), CED = CED, 
      aa = aa, CES = CES, ces.ans)
  bb <- root.out$root
  if (bb < 0) {
  }
  return(bb)
}


#' Determine values for the CED; categorical model
#' @param ans.all ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with matrix of CED, response matrix, pi background matrix and gr.txt
#' @export
f.ced.cat <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.ced.cat")
  
  ans.all <- f.basics.cat(ans.all, track = track)
  
  with(ans.all, {
        
        CED.matr <- matrix(nrow = nr.gr, ncol = nth)
        special <- FALSE
        
        if ((dtype %in% c(4, 6)) && (model.ans %in% c(12:15, 22:25))) {
          if (nr.aa > 1 && nr.bb == 1) {
            CED.matr <- matrix(regr.par[2:(nr.gr + 1)], ncol = 1)
            special <- TRUE
          }
        }
        
        if (ces.ans %in% c(1:3) & model.type == 2 & model.ans %in% c(12:15, 22:25)) {
          
          if (model.ans < 16) 
            model.ans <- model.ans - 10
          if (model.ans > 21) 
            model.ans <- model.ans - 5
          ans.all$regr.par <- f.bb.cat(ans.all)
          ans.all$model.ans <- model.ans
          regr.par.matr <- f.pars(ans.all)$regr.par.matr
        }
        
        
        if (model.type == 1)
          for (gr in (1:nr.gr))
            CED.matr[gr, 1] <- f.inv.bin(model.ans = model.ans, 
                regr.par = regr.par.matr[gr, ], 
                CES = CES.all[gr], ces.ans = ces.ans, track = track)
        
        
        if (model.type == 2) {
          for (gr in (1:nr.gr)) {
            if (!special) 
              CED.matr[gr, ] <- f.inv.con(model.ans, 
                  params = regr.par.matr[gr, ],
                  CES = CES.con.matr[gr, ], dtype = dtype)
          }
        }
        
        if (ces.ans == 1) {
          
          pi.bg.up <- matrix(0, nr.gr, (nth + 1))
          for (k in 1:nth) {
            dum <- nth + 1 - k
            pi.bg.up[, dum] <- pi.bg.matr[, dum] + pi.bg.up[, (dum + 1)]
          }
          pi.bg.up <- pi.bg.up[, 1:nth]
          if (!special) 
            CED.matr[pi.bg.up > 0.5] <- -Inf
          model.lst <- c(4, 5, 9, 10, 14, 15)
          th.cumu <- exp(cumsum(th.par))
          if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 24, 25) & model.type == 2) 
            for (jj in (1:nth)) 
              for (gr in (1:nr.gr)) {
                if (regr.par.matr[gr, 1] * regr.par.matr[gr, 3] > th.cumu[jj]) 
                  if (!special) 
                    CED.matr[gr, jj] <- Inf
              }
        }
        
        CED.lst <- list(CED.matr = CED.matr, response.matr = response.matr, 
            pi.bg.matr = pi.bg.matr, gr.txt = ans.all$gr.txt)
        
        if (track)	print("f.ced.cat : END")
        
        return(CED.lst)
        
      })
  
}

#' compute basic statistics of a categorical model
#' @param ans.all ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return update ans.all
#' @export
f.basics.cat <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.basics.cat")
  
  with(ans.all, {
        
        nrp <- length(regr.par) - 1
        nth <- ifelse(model.type > 1, length(th.par), 1)
        response.matr <- matrix(NA, nrow = nr.gr, ncol = nth)
        pi.bg.matr <- matrix(NA, nrow = nr.gr, ncol = nth)
        CES.con.matr <- matrix(NA, nrow = nr.gr, ncol = nth)
        CES.all <- rep(CES, nr.gr)
        ans.all <- f.pars(ans.all)
        regr.par.matr <- ans.all$regr.par.matr
        
        if (ces.ans > 1) 
          
          if (!decr.zz) {
            if(track)	cat("\n\nATTENTION:  quantal dose-response is treated as decreasing\n\n")
            CES.all <- -CES.all
          }
        
        switch(as.character(model.type),
            
            '1' = {
              
              for (gr in (1:nr.gr)) 
                pi.bg.matr[gr, 1] <- f.expect.bin(model.ans = model.ans, 
                    x = 0, regr.par = regr.par.matr[gr, ], CES = CES, 
                    ces.ans = ces.ans, track = track)
              
              switch(as.character(ces.ans),
                  '1' = response.matr[, ] <- 0.5,
                  '2' = response.matr <- pi.bg.matr + CES.all,
                  '3' = response.matr <- pi.bg.matr + CES.all * (1 - pi.bg.matr),
                  '4' = response.matr <- pnorm(CES.all + qnorm(pi.bg.matr)),
                  '5' = response.matr <- pi.bg.matr * (CES.all + 1)
              )
            },
            
            '2' = {
              
              if (model.ans %in% c(12:15, 22:25) && ces.ans < 4) {
                # in ans.all.tmp, change model.ans from
                # 12:15 -> 2:5, 22:25 -> 17:20 
                ans.all.tmp <- f.model.bb(ans.all, track = track)
                model.ans <- ans.all.tmp$model.ans
                regr.par.matr <- f.pars(ans.all.tmp)$regr.par.matr
              }
              
              for (gr in (1:nr.gr)) {
                
                uu.bg <- regr.par.matr[gr, 1]
                zz.bg <- log(uu.bg)
                pi.bg.matr[gr, ] <- f.expect.cat(model.type = 2, 
                    model.ans, x = 0, regr.par = regr.par.matr[gr, ], th.par, 
                    sig.par, CES = CES, ces.ans = ces.ans, dtype = dtype, 
                    twice = twice, track = track)
                
                if (ces.ans == 1) {
                  zz.mu <- cumsum(th.par)
                  if (decr.zz) 
                    zz.mu[zz.mu > logb(regr.par.matr[gr, 1])] <- NA
                  if (!decr.zz) 
                    zz.mu[zz.mu < logb(regr.par.matr[gr, 1])] <- NA
                  response.matr[, ] <- 0.5
                }
                
                if (ces.ans > 1) {
                  if (ces.ans == 2) 
                    response.matr[gr, ] <- pi.bg.matr[gr, ] + CES.all[gr]
                  if (ces.ans == 3) 
                    response.matr[gr, ] <- pi.bg.matr[gr, ] + CES.all[gr] * (1 - pi.bg.matr[gr, ])
                  snd <- qnorm(response.matr[gr, ])
                  zz.mu <- cumsum(th.par) - (snd * sig.par)
                }
                
                CES.con.matr[gr, ] <- 1 - exp(zz.mu)/uu.bg
                if (!decr.zz)	CES.con.matr[, gr] <- -CES.con.matr[, gr]
              }
              
            }
        )
        
        ans.all$response.matr <- response.matr
        ans.all$pi.bg.matr <- pi.bg.matr
        ans.all$CES.all <- CES.all
        ans.all$CES.con.matr <- CES.con.matr
        
        if (track)	print("f.basics.cat : END")
        
        return(ans.all)
      })
  
}

#' computed expected value for categorical model
#' @param model.type model type
#' @param model.ans type or response fitted
#' @param x observations
#' @param regr.par parameters
#' @param th.par parameter theta
#' @param sig.par parameter sigma (standard deviations)?
#' @param fct1 numeric, value for first parameter
#' @param fct2 numeric, value for second parameter
#' @param CES CES value
#' @param CES.cat CES cat
#' @param latent latent model
#' @param ttt ttt
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param fct3 numeric, value for third parameter
#' @param ces.ans type of benchmark response
#' @param decr.zz is decreasing
#' @param dtype response data type
#' @param fct4 numeric, value for fourth parameter
#' @param fct5 numeric, value for fifth parameter
#' @param fct1.txt text for parameter 1
#' @param fct2.txt text for parameter 2
#' @param fct3.txt text for parameter 3
#' @param fct4.txt text for parameter 4
#' @param fct5.txt text for parameter 5
#' @param covar.txt text for covar
#' @param fct3.ref parameter 3 of reference
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return expected value
#' @export
f.expect.cat <- function (model.type = 2, model.ans, x, regr.par = 0, th.par = 0, 
    sig.par = 0, fct1 = 1, fct2 = 1, CES = 0, CES.cat = 1, 
    latent = F, ttt = 0, twice, fct3 = 1, ces.ans, decr.zz = T, 
    dtype, fct4 = 1, fct5 = 1, fct1.txt = "", fct2.txt = "", 
    fct3.txt = "", fct4.txt = "", fct5.txt = "", covar.txt = "", 
    fct3.ref, track = FALSE) {
  
  if (track)	print("f.expect.cat")
  
  nth <- length(th.par)
  ans.all.tmp <- list()
  ans.all.tmp$model.ans <- model.ans
  ans.all.tmp$model.type <- model.type
  ans.all.tmp$dtype <- dtype
  ans.all.tmp$quick.ans <- 1
  ans.all.tmp$regr.par <- regr.par
  ans.all.tmp$th.par <- th.par
  ans.all.tmp$nth < nth
  ans.all.tmp$sig.par <- sig.par
  ans.all.tmp$twice <- twice
  ans.all.tmp$fct1 <- fct1
  ans.all.tmp$fct2 <- fct2
  ans.all.tmp$fct3 <- fct3
  ans.all.tmp$fct4 <- fct4
  ans.all.tmp$fct5 <- fct5
  ans.all.tmp$fct1.txt <- fct1.txt
  ans.all.tmp$fct2.txt <- fct2.txt
  ans.all.tmp$fct3.txt <- fct3.txt
  ans.all.tmp$fct4.txt <- fct4.txt
  ans.all.tmp$fct5.txt <- fct5.txt
  ans.all.tmp$covar.txt <- covar.txt
  ans.all.tmp$CES <- CES
  ans.all.tmp$CES.cat <- CES.cat
  ans.all.tmp$ces.ans <- ces.ans
  ans.all.tmp$decr.zz <- decr.zz
  increase <- 1
  
  if (decr.zz)	increase <- -1
  
  if (model.ans %in% c(12:15, 22:25)) {
    ans.all.tmp$increase <- increase
    ans.all.tmp$nr.var <- max(fct3)
    ans.all.tmp$nr.cc <- max(fct4)
    ans.all.tmp$nr.dd <- max(fct5)
    ans.all.tmp <- f.model.bb(ans.all.tmp)
    regr.par <- ans.all.tmp$regr.par
    model.ans <- ans.all.tmp$model.ans
  }
  
  uu <- f.expect.con(model.ans, x, regr.par, fct1, fct2, fct5 = fct5, 
      ttt = ttt, twice = twice, CES = -abs(CES), increase = increase, track = track)
  
  zz <- logb(uu)
  
  if (length(x) < 200) 
    if (latent == T) {
      return(zz)
    }
  
  th.lvm <- NA
  if (max(fct3) > 1) {
    th.lvm <- f.th.lvm(th.par, fct3, fct3.ref, track = track)
    pipi <- matrix(, length(x), 1)
    nr.th <- length(th.lvm)
    fct3 <- as.numeric(factor(fct3))
    th0 <- rep(0, length(x))
    for (ii in 1:nr.th) {
      th0 <- th0 + th.lvm[ii] * (fct3 == ii)
    }
    pipi[, 1] <- pnorm((th0 - zz)/sig.par)
  }else {
    nth <- length(th.par)
    th <- cumsum(th.par)
    pipi <- matrix(, length(x), nth)
    right <- 0
    for (k in 1:nth) {
      dum <- nth - k + 1
      pipi[, dum] <- pnorm((th[dum] - zz)/sig.par) - right
      right <- right + pipi[, dum]
    }
  }
  
  if (track)	print("f.expect.cat: END")
  
  return(pipi)
}

#' change model.ans
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return updated ans.all
#' @export
f.model.bb <- function (ans.all, track = FALSE) 
{
  
  if (track)	print("f.model.bb")
  
  with(ans.all, {
        
        if (model.ans %in% c(12:15, 22:25)) {
          ans.all$regr.par <- f.bb.cat(ans.all = ans.all, track = track)
          if (model.ans %in% 12:15) 
            ans.all$model.ans <- model.ans - 10
          if (model.ans %in% 22:25) 
            ans.all$model.ans <- model.ans - 5
          ans.all$CES <- -CES
        }
        
        if (track)	print("f.model.bb: END")
        return(ans.all)
        
      })
  
}

#' ?
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return ?
#' @export
f.bb.cat <- function (ans.all, track = FALSE) 
{
  if (track)	print("f.bb.cat")
  
  with(ans.all, {
        
        if (ces.ans == 1) 
          return(f.par.u(model.ans = model.ans, regr.par = regr.par, th.par = th.par, 
                  nr.aa = max(fct1), nr.bb = max(fct2), CES.cat = CES.cat, track = track))
        nr.aa <- max(fct1)
        nth <- length(th.par)
        regr.par.matr <- f.pars(ans.all)$regr.par.matr
        nr.gr <- length(regr.par.matr[, 1])
        response.matr <- matrix(NA, nrow = nr.gr, ncol = nth)
        pi.bg.matr <- matrix(nrow = nr.gr, ncol = nth)
        
        # function 1
        f.uniroot.aa <- function(bb, par.rest, th.par, model.ans, CED, CES) {
          
          if (track)	print("f.uniroot.aa")
          
          f.aa <- function(aa, bb, par.rest, th.par, model.ans, CED, CES) {
            zz.bg <- log(aa)
            pi.bg <- pnorm(cumsum(th.par) - zz.bg, 0, sig.par)
            resp <- pi.bg + CES * (1 - pi.bg)
            zz.mu <- cumsum(th.par) - qnorm(resp, 0, sig.par)
            regr.par.tmp <- c(aa, bb, par.rest)
            uu.mu <- f.expect.con(model.ans, CED, regr.par.tmp, increase = increase, track = track)
            return(zz.mu - log(uu.mu))
          }
          
          try1 <- 0
          try2 <- 0
          aa.low <- 0.01
          aa.upp <- 10
          go.on <- T
          count <- 1
          while (go.on) {
            try1 <- f.aa(aa.low, bb, par.rest, th.par, model.ans, CED, CES)
            try2 <- f.aa(aa.upp, bb, par.rest, th.par, model.ans, CED, CES)
            if (is.na(try1)) 
              return(NA)
            if (is.na(try2)) 
              return(NA)
            if (sign(try1) == sign(try2)) 
              go.on <- T
            else go.on <- F
            count <- count + 1
            if (count > 5)	go.on <- F
          }
          if (count <= 5) {
            root.out <- uniroot(f.aa, interval = c(aa.low, aa.upp), 
                bb, par.rest, th.par, model.ans, CED, CES)
            return(root.out$root)
          }else return(1)
        }
        
        for (gr in 1:nr.gr) {
          
          uu.bg <- regr.par.matr[gr, 1]
          zz.bg <- log(uu.bg)
          pi.bg.matr[gr, ] <- pnorm(cumsum(th.par) - zz.bg, 0, sig.par)
          
          if (ces.ans == 2) {
            response.matr[gr, ] <- pi.bg.matr[gr, ] + CES
            zz.mu <- cumsum(th.par) - qnorm(response.matr[gr], 0, sig.par)
          }
          
          if (ces.ans == 3) {
            response.matr[gr, ] <- pi.bg.matr[gr, ] + CES * (1 - pi.bg.matr[gr, ])
            zz.mu <- cumsum(th.par) - qnorm(response.matr[gr], 0, sig.par)
          }
          
          CES.con <- exp(zz.mu)/exp(zz.bg) - 1
          if (!decr.zz)	CES.con <- -CES.con
          CED <- regr.par.matr[gr, 2]
          cc <- NA
          dd <- NA
          nrp <- length(regr.par.matr[1, ])
          if (model.ans %in% c(13, 23)) 
            dd <- regr.par.matr[gr, nrp]
          if (model.ans %in% c(14, 24)) 
            cc <- regr.par.matr[gr, nrp]
          if (model.ans %in% c(15, 25)) 
            cc <- regr.par.matr[gr, nrp - 1]
          if (model.ans %in% c(15, 25)) 
            dd <- regr.par.matr[gr, nrp]
          if (is.na(CES.con)) {
            if(track)	cat("\nCES.con in f.bb.cat is NA\n")
            bb <- NA
          }else bb <- f.bb.con(model.ans, cc = cc, dd = dd, CED, CES.con)
          if (gr == 1) 
            bb.ref <- bb
          if (max(fct1) > 1 & max(fct2) == 1) {
            if (gr == 1) {
              regr.par[nr.aa + 1] <- bb
            }else {
              if (model.ans %in% 12:15) 
                model.ans.tmp <- model.ans - 10
              if (model.ans %in% 22:25) 
                model.ans.tmp <- model.ans - 5
              par.rest <- regr.par.matr[gr, -(1:2)]
              aa <- f.uniroot.aa(bb.ref, par.rest, th.par, model.ans.tmp, CED, CES)
              regr.par[gr] <- aa
            }
          }
          else {
            regr.par[nr.aa + gr] <- bb
          }
        }
        
        if (track)	print("bb.cat: END")
        
        return(regr.par)
        
      })
  
}

#' ? 
#' 
#' called in the f.bb.cat function
#' @param model.ans type of response
#' @param regr.par regression parameter
#' @param th.par theta parameter
#' @param nr.aa number of a
#' @param nr.bb number of b
#' @param CES.cat ?
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return regr.par.u
#' @export
f.par.u <- function (model.ans, regr.par = 0, th.par = 0, nr.aa = 1, nr.bb = 1, CES.cat, track = FALSE) {
  
  if (track)	print("f.par.u")
  
  if (model.ans %in% 12:15) 
    model.ans.tmp <- model.ans - 10
  else if (model.ans %in% 22:25) 
    model.ans.tmp <- model.ans - 15
  else return("f.par.u has been called for inadequate model")
  
  nr.CED <- max(nr.aa, nr.bb)
  th <- cumsum(th.par)
  if (mode(CES.cat) == "NULL" & track) 
    cat("\nATTENTION: CES.cat unknown, retry by choosing model again\n")
  th.CES <- th[CES.cat]
  par3 <- regr.par[nr.aa + nr.CED + 1]
  par4 <- regr.par[nr.aa + nr.CED + 2]
  regr.par.u <- 1
  
  if (nr.aa == 1 & nr.bb == 1) {
    
    regr.par.u[1] <- regr.par[1]
    th.CES <- th.CES - logb(regr.par[1])
    
    switch(as.character(model.ans.tmp), 
        
        '2' = regr.par.u[2] <- th.CES/regr.par[2],
        
        '3' = {
          regr.par.u[2] <- th.CES/(regr.par[2]^par3)
          regr.par.u[3] <- par3
        },
        
        '4' = {
          regr.par.u[2] <- -(1/regr.par[2]) * logb((exp(th.CES) - par3)/(1 - par3))
          regr.par.u[3] <- par3
        }, 
        
        '5' = {
          regr.par.u[2] <- -(1/(regr.par[2]^par4)) * logb((exp(th.CES) - par3)/(1 - par3))
          regr.par.u[3] <- par3
          regr.par.u[4] <- par4
        }, 
        
        '7' = {
          tmp <- exp(th.CES)
          CED.1 <- regr.par[2]
          regr.par.u[2] <- CED.1 * tmp/(1 - tmp)
        }, 
        
        '8' = {
          tmp <- exp(th.CES)
          CED.1 <- regr.par[2]
          regr.par.u[2] <- CED.1 * (tmp/(1 - tmp))^(1/par3)
          regr.par.u[3] <- par3
        }, 
        
        '9' = {
          CED.1 <- regr.par[2]
          tmp <- exp(th.CES)
          bb <- CED.1 * ((par3 - tmp)/(tmp - 1))
          regr.par.u[2] <- bb
          regr.par.u[3] <- par3
        }, 
        
        '10' = {
          CED.1 <- regr.par[2]
          tmp <- exp(th.CES)
          bb <- CED.1 * ((par3 - tmp)/(tmp - 1))^(1/par4)
          regr.par.u[2] <- bb
          regr.par.u[3] <- par3
          regr.par.u[4] <- par4
        }
    )
  }
  
  if (nr.aa > 1 & nr.bb == 1) {
    
    regr.par.u[1] <- regr.par[1]
    th.tmp <- (th.CES - logb(regr.par[1]))
    par3 <- regr.par[nr.aa + nr.CED]
    par4 <- regr.par[nr.aa + nr.CED + 1]
    
    switch(as.character(model.ans.tmp),
        
        '2' = {
          CED.1 <- regr.par[2]
          bb <- th.tmp/CED.1
          for (ii in 2:(nr.aa))
            regr.par.u[ii] <- exp(th.CES - bb * regr.par[ii + 1])
          regr.par.u[nr.aa + 1] <- bb	
        }, 
        
        '3' = {
          CED.1 <- regr.par[2]
          bb <- th.tmp/(CED.1^par3)
          for (ii in 2:(nr.aa)) 
            regr.par.u[ii] <- exp(th.CES - bb * (regr.par[ii + 1]^par3))
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
        }, 
        
        '4' = {
          CED.1 <- regr.par[2]
          bb <- -log((exp(th.tmp) - par3)/(1 - par3))/CED.1
          for (ii in 2:nr.aa)
            regr.par.u[ii] <- exp(th.CES - logb(par3 - (par3 - 1) * exp(-bb * regr.par[ii + 1])))
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
        }, 
        
        '5' = {
          CED.1 <- regr.par[2]
          bb <- -log((exp(th.tmp) - par3)/(1 - par3))/(CED.1^par4)
          for (ii in 2:nr.aa)
            regr.par.u[ii] <- exp(th.CES - logb(par3 - (par3 - 1) * exp(-bb * regr.par[ii + 1]^par4)))
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
          regr.par.u[nr.aa + 3] <- par4
        }, 
        
        '7' = {
          CED.1 <- regr.par[2]
          th.tmp <- exp(th.tmp)
          bb <- CED.1 * th.tmp/(1 - th.tmp)
          for (ii in 2:nr.aa) {
            CED.dum <- regr.par[ii + 1]
            dum <- CED.dum/(bb + CED.dum)
            regr.par.u[ii] <- exp(th.CES - log(1 - dum))
          }
          regr.par.u[nr.aa + 1] <- bb
        }, 
        
        '8' = {
          CED.1 <- regr.par[2]
          th.tmp <- exp(th.tmp)
          bb <- CED.1 * (th.tmp/(1 - th.tmp))^(1/par3)
          for (ii in 2:nr.aa) {
            CED.dum <- regr.par[ii + 1]
            dum <- CED.dum^par3/(bb^par3 + CED.dum^par3)
            regr.par.u[ii] <- exp(th.CES - log(1 - dum))
          }
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
        }, 
        
        '9' = {
          CED.1 <- regr.par[2]
          th.tmp <- exp(th.tmp)
          bb <- CED.1 * ((par3 - th.tmp)/(th.tmp - 1))
          for (ii in 2:nr.aa) {
            CED.dum <- regr.par[ii + 1]
            dum <- CED.dum/(bb + CED.dum)
            regr.par.u[ii] <- exp(th.CES - log(1 + (par3 - 1) * dum))
          }
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
        }, 
        
        '10' = {
          CED.1 <- regr.par[2]
          CED.1 <- regr.par[2]
          th.tmp <- exp(th.tmp)
          bb <- CED.1 * (((par3 - th.tmp)/(th.tmp - 1)))^(1/par4)
          for (ii in 2:nr.aa) {
            CED.dum <- regr.par[ii + 1]
            dum <- CED.dum^par4/(bb^par4 + CED.dum^par4)
            regr.par.u[ii] <- exp(th.CES - log(1 + (par3 - 1) * dum))
          }
          regr.par.u[nr.aa + 1] <- bb
          regr.par.u[nr.aa + 2] <- par3
          regr.par.u[nr.aa + 3] <- par4
        }
    )
  }
  
  if (nr.bb > 1) {
    
    regr.par.u <- regr.par
    switch(as.character(model.ans.tmp), 
        
        '2' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1)
              aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            CED.tmp <- regr.par[nr.aa + jj]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- log(tmp)/CED.tmp
          }
        }, 
        
        '3' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            CED.tmp <- regr.par[nr.aa + jj]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- log(tmp)/CED.tmp^par3
          }
        }, 
        
        '4' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            CED.tmp <- regr.par[nr.aa + jj]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- -(1/CED.tmp) * logb((exp(th.CES - log(aa.tmp)) - par3)/(1 - par3))
          }
        }, 
        
        '5' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            CED.tmp <- regr.par[nr.aa + jj]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- -(1/CED.tmp^par4) * 
                logb((exp(th.CES - log(aa.tmp)) - par3)/(1 - par3))
          }
        }, 
        
        '7' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            CED.tmp <- regr.par[nr.aa + jj]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- regr.par[nr.aa + jj] * tmp/(1 - tmp)
          }
        }, 
        
        '8' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- regr.par[nr.aa + jj] * (tmp/(1 - tmp))^(1/par3)
          }
        }, 
        
        '9' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- regr.par[nr.aa + jj] * ((par3 - tmp)/(tmp - 1))
          }
        }, 
        
        '10' = {
          for (jj in 1:nr.bb) {
            if (nr.aa > 1) aa.tmp <- regr.par[jj] else aa.tmp <- regr.par[1]
            tmp <- exp(th.CES - log(aa.tmp))
            regr.par.u[nr.aa + jj] <- regr.par[nr.aa + jj] * (((par3 - tmp)/(tmp - 1)))^(1/par4)
          }
        }
    )
  }
  
  if(track)	print("f.par.u : END")
  
  return(regr.par.u)
  
}


#' Define initial parameter values for the chosen model, categorical response
#' 
#' from f.start.cat of the proast61.3 package, only: tmp.quick == FALSE
#' @param ans.all list, with all results that were obtained during the analysis
#' @param adjust boolean, if TRUE the start values of the model parameters will
#' be adjusted; default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.start.cat <- function (ans.all, adjust = FALSE, track = FALSE) {
  
  if (track)	print("f.start.cat")
  
  ans.all$adjust <- adjust
  ans.all$tmp.quick <- FALSE
  
  if (ans.all$decr.zz) 
    ans.all$increase <- -1
  
  with(ans.all, {
        sig.start <- 1
        th.0.start <- 0
        ans.all$th.0.start <- th.0.start
        ans.all$sig.start <- sig.start
        th.start <- c(th.0.start, rep(-sig.start, (nth - 1)))
        if (dtype == 2) 
          th.start <- th.0.start
        if (model.ans %in% c(12:15, 22:25)) 
          if (nr.aa > 1 && nr.bb > 1 && fct1.no != fct2.no) {
            cat("\nLVM models in terms of CED are not implemented for different covariates on a and b\n")
            cat("adjust covariate before proceeding\n")
            return(ans.all)
          }
        
        if (dtype %in% 2:3) {
          lmp.lst <- f.lump.cat(x, y, nn, fct1, fct2, dtype, 
              twice)
          xlmp <- lmp.lst$xlmp
          ylmp <- lmp.lst$ylmp
          nlmp <- lmp.lst$nlmp
          fct1.lmp <- lmp.lst$fct1.lmp
          fct2.lmp <- lmp.lst$fct2.lmp
          ylmp <- ylmp + 0.01 * (ylmp == 0) - 0.01 * (ylmp == 
                1)
        } else {
          xlmp <- x
          ylmp <- y
          nlmp <- nn
          fct1.lmp <- fct1
          fct2.lmp <- fct2
        }
        
        if (dtype %in% c(4, 6, 84)) {
          ylmp.corr <- ylmp + 0.01 * (ylmp == 0) - 0.01 * (ylmp == 
                1)
          dt.fr <- data.frame(yyy = logb(1 - ylmp.corr), dose.tmp = xlmp)
          res.lm <- lm(yyy ~ dose.tmp, dt.fr)
          bb <- 1/res.lm[[1]][2]
          if (!tmp.quick && !is.finite(bb)) {
            cat("\nATTENTION: no dose-response could be detected, check your data\n")
            bb <- 1
          }
        }
        
        if(!adjust){
          
          ans.all.tmp <- ans.all
          ans.all.tmp$x <- xlmp
          ans.all.tmp$yy <- exp(ylmp)
          ans.all.tmp <- f.start.con(ans.all.tmp, track = track)
          regr.start <- ans.all.tmp$regr.par
          aa <- regr.start[1:nr.aa]
          aa[aa < 0.1] <- 1
          aa[aa < exp(th.0.start)] <- exp(th.0.start) * 1.1
          regr.start[1:nr.aa] <- aa
          ans.all$regr.start <- regr.start
          ans.all$regr.par <- regr.start
          bb <- regr.start[(nr.aa + 1):(nr.aa + nr.bb)]
          nr.CED <- max(nr.aa, nr.bb)
          if (model.ans %in% 12:15) {
            CED <- mean(x)/2
            if (dtype == 3) 
              CED <- mean(x)
            if (ces.ans %in% 2:3) 
              CED <- 10 * CES * CED
          }
          if (model.ans %in% 22:25) {
            if (ces.ans == 1) 
              CED <- mean(x)
            if (ces.ans > 1) 
              CED <- mean(x)/7
            if (ces.ans %in% 2:3) 
              CED <- 10 * CES * CED
          }
          
          switch(model.ans, {
                if (nr.bb > 1) {
                  cat("\n\nthis model can not be fitted for factor-dependent parameter b\n")
                  return(invisible())
                }
              }, {
                regr.start[(nr.aa + 1):(nr.aa + nr.bb)] <- -abs(bb)
              }, {
                regr.start[(nr.aa + 1):(nr.aa + nr.bb)] <- -abs(bb)
              }, {
                cc <- regr.start[nr.aa + nr.bb + 1]
                if (cc > 1) cc <- 1/cc
                regr.start[nr.aa + nr.bb + 1] <- cc
              }, {
                cc <- regr.start[nr.aa + nr.bb + 1]
                if (cc > 1) cc <- 1/cc
                regr.start[nr.aa + nr.bb + 1] <- cc
              }, {
                cc <- 0.01 * ytmp[top]/regr.start[1]
                if (!is.finite(bb)) bb <- xtmp[5]/100
              }, {
                regr.start[(nr.aa + 1):(nr.aa + nr.bb)] <- -abs(bb)
              }, {
                regr.start[(nr.aa + 1):(nr.aa + nr.bb)] <- -abs(bb)
              }, {
                cc <- regr.start[nr.aa + nr.bb + 1]
                if (cc > 1) cc <- 1/cc
                regr.start[nr.aa + nr.bb + 1] <- cc
              }, {
                cc <- regr.start[nr.aa + nr.bb + 1]
                if (cc > 1) cc <- 1/cc
                regr.start[nr.aa + nr.bb + 1] <- cc
              }, {
                eps <- 1e-06
                if (dtype == 2) {
                  regr.start <- as.numeric(y)
                  regr.start <- regr.start + eps * (regr.start == 
                        0) - eps * (regr.start == 1)
                  cat("\nNOTE:  Saturated model may not be applicable for dtype = 2\n\n")
                }
                if (dtype == 6) {
                  kk.tmp <- tapply(kk, x, sum)
                  nn.tmp <- tapply(nn, x, sum)
                  x.tmp <- as.numeric(tapply(x, x, mean))
                  regr.start <- as.numeric(kk.tmp/nn.tmp)
                  regr.start <- regr.start + eps * (regr.start == 
                        0) - eps * (regr.start == 1)
                  if (max(nr.aa, nr.bb) > 1) cat("\nATTENTION: covariates not yet implemented for nested quantal data\n")
                }
                if (dtype == 4) {
                  regr.start <- as.numeric(kk/nn)
                  regr.start <- regr.start + eps * (regr.start == 
                        0) - eps * (regr.start == 1)
                }
              }, {
                regr.start <- c(regr.start[1], rep(CED, nr.CED))
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED))
              }, {
                dd <- 0.5
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    dd)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), dd)
              }, {
                th.cum <- cumsum(th.start)
                th.CES <- th.cum[CES.cat]
                cc <- exp(th.CES)/10
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    cc)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), cc)
              }, {
                th.cum <- cumsum(th.start)
                th.CES <- th.cum[CES.cat]
                cc <- exp(th.CES)/10
                dd <- 0.5
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    cc, dd)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), cc, dd)
              }, cat(""), {
                bb <- mean(x)/3
                if (dtype == 3) bb <- mean(x)
                regr.start <- c(aa, rep(bb, nr.bb))
              }, {
                bb <- mean(x)/3
                if (dtype == 3) bb <- mean(x)
                dd <- 0.5
                regr.start <- c(aa, rep(bb, nr.bb), dd)
              }, {
                bb <- mean(x)/3
                if (dtype == 3) bb <- mean(x)
                cc <- 0.1
                regr.start <- c(aa, rep(bb, nr.bb), cc)
              }, {
                bb <- mean(x)/3
                if (dtype == 3) bb <- mean(x)
                cc <- 0.1
                dd <- 0.5
                regr.start <- c(aa, rep(bb, nr.bb), cc, dd)
              }, cat("model 21"), {
                regr.start <- c(regr.start[1], rep(CED, nr.CED))
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED))
              }, {
                dd <- 0.5
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    dd)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), dd)
              }, {
                th.cum <- cumsum(th.start)
                th.CES <- th.cum[CES.cat]
                cc <- exp(th.CES)/10
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    cc)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), cc)
              }, {
                dd <- 0.5
                th.cum <- cumsum(th.start)
                th.CES <- th.cum[CES.cat]
                cc <- exp(th.CES)/10
                regr.start <- c(regr.start[1], rep(CED, nr.CED), 
                    cc, dd)
                if (nr.aa > 1 && nr.bb > 1) regr.start <- c(regr.start[1:nr.aa], 
                      rep(CED, nr.CED), cc, dd)
              })
          
          par.start <- c(regr.start, th.start, sig.start)
          nrp <- length(regr.start)
          npar <- length(par.start)
          ans.all$par.start <- par.start
          ans.all$regr.start <- regr.start
          ans.all$th.start <- th.start
          ans.all$sig.start <- sig.start
          ans.all$npar <- npar
          ans.all$nrp <- nrp
          ans.all$nth <- nth
          ans.all <- f.constr.con(ans.all)
          if (dtype == 6) {
            par.start <- c(alfa.start, par.start)
            if (!(model.ans == 14 & model.type == 1)) 
              par.start[1] <- alfa.start
            npar <- length(par.start)
            ans.all$par.start <- par.start
            ans.all$npar <- npar
          }
          loglik.old <- NA
          
          loglik.old <- -f.lik.cat(ans.all$par.start, x, 
              y, kk, nn, dtype, fct1, fct2, ans.all$nrp, 
              ans.all$nth, nr.aa, nr.bb, model.ans, model.type, 
              ttt = ttt, twice = twice, fct3 = fct3, 
              ces.ans = ces.ans, CES = CES, CES.cat = CES.cat, 
              decr.zz = decr.zz, alfa.length = alfa.length, 
              fct3.ref = fct3.ref, kk.tot = kk.tot, nn.tot = nn.tot)
          if(track)
            cat("\n\n log-likelihood value: ", loglik.old, "\n\n")
          
          
          ans.all$loglik.old <- loglik.old
          
        } else { # if adjust == TRUE
          
          loglik.old <- 0
          npar <- length(ans.all$par.start)
          
          ans.all$par.start <- startValues
          par.lst <- f.split.par(ans.all$par.start, nrp, nth, 
              dtype)
          regr.start <- par.lst$regr.par
          th.start <- par.lst$th.par
          sig.start <- par.lst$sig.par
          
          loglik.new <- NA
          
          loglik.new <- -f.lik.cat(ans.all$par.start, 
              x, y, kk, nn, dtype, fct1, fct2, nrp, nth, 
              nr.aa, nr.bb, model.ans, model.type, ttt = ttt, 
              twice = twice, fct3 = fct3, 
              ces.ans = ces.ans, CES = CES, CES.cat = CES.cat, 
              decr.zz = decr.zz, alfa.length = alfa.length, 
              fct3.ref = fct3.ref, kk.tot = kk.tot, nn.tot = nn.tot)
          
          if(is.na(loglik.new))
            ans.all$errorAdjustStartValues <- TRUE
          
          
          if (model.type == 2) 
            if (model.ans %in% c(14:15, 24:25)) {
              nr.CED <- max(nr.aa, nr.bb)
              cc <- regr.start[nr.aa + nr.bb + 1]
              th.cum <- cumsum(th.start)
              th.CES <- th.cum[CES.cat]
              
              # MV added this line
              if(cc > 0)
                if (CES == 0 & th.CES < logb(cc)) {
                  warning(paste("\n\nATTENTION: parameter c should be smaller than theta at CED\n\n                          Make c smaller than ", 
                          exp(th.CES), "\n\n"))
                  ans.tmp <- 0
                }
            }
          
          loglik.old <- loglik.new
          
        }
        
        ans.all$regr.start <- regr.start
        
        if (twice) 
          max.lev <- max(max(fct1), max(fct2))
        if (!twice) 
          max.lev <- max(fct1) * max(fct2)
        if (dtype == 4) 
          max.lev <- max.lev * nth
        if (dtype == 84) 
          max.lev <- max(as.numeric(factor(ttt)))
        ans.all$max.lev <- max.lev
        
        if (track) print("f.start.cat : END")
        
        return(ans.all)
        
      })
}



#' Compare model fitted with and without alfa
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all 
#' @export
f.dtype6.mn <- function (ans.all, track = FALSE) {
  
  if (track) print("f.dtype6.mn")
  
  ans.all.tmp <- ans.all
  ans.all.tmp$plot.type <- 0
  ans.all.tmp$model.ans <- 14
  ans.all.tmp$model.type <- 1
  ans.all.tmp$model.name <- "full model"
  ans.all.tmp <- f.start.bin(ans.all = ans.all.tmp, track = track)
  ans.all.tmp$lb[1] <- 1e+05
  ans.all.tmp$ub[1] <- 1e+05
  if(track)
    cat("\n fitting full model, without alfa .....\n")
  ans.all.tmp <- f.nlminb(ans.all.tmp, track)
  loglik.1 <- ans.all.tmp$loglik
  ans.all.tmp$lb[1] <- 1e-06
  ans.all.tmp$ub[1] <- Inf
  ans.all.tmp <- f.start.bin(ans.all.tmp, track = track)
  if(track)
    cat("\n fitting full model, with alfa (find group means) .....\n")
  ans.all.tmp <- f.nlminb(ans.all.tmp, track)
  loglik.2 <- ans.all.tmp$loglik
  Pvalue <- f.P(loglik.1, loglik.2, 1, track = track)
  if(track)
    cat("\nThe P-value for adding parameter alfa to the full model is:", 
        Pvalue, "\n\n")
  ans.all.tmp$fitted <- F
  ans.all.tmp$pi.full <- ans.all.tmp$MLE[-(1:ans.all$alfa.length)]
  ans.all.tmp$alfa.start <- ans.all.tmp$MLE[1]
  ans.all.tmp$plot.type <- ans.all$plot.type
  ans.all.tmp$model.ans <- ans.all$model.ans
  ans.all.tmp$model.type <- ans.all$model.type
  ans.all.tmp$model.name <- ans.all$model.name
  ans.all.tmp$Pvalue.alfa <- Pvalue
  
  if (track) print("f.dtype6.mn:  END")
  
  return(ans.all.tmp)
}



#' Calculate p-value for likelihood test 
#' @param loglik.1 numeric, first likelihood value 
#' @param loglik.2 numeric, second likelihood value
#' @param df integer, degrees of freedom
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric, p-value for likelihood test 
#' @export
f.P <- function (loglik.1, loglik.2, df, track = FALSE) {
  
  if (df <= 0) {
    cat("\nf.P:  nonpositive df !")
    cat("\ndf=", df)
    cat("\nloglik 1=", loglik.1)
    cat("\nloglik 2=", loglik.2)
    cat("\n")
  }
  
  if (df == 0) 
    df <- 1
  
  if (track) 
    print(c(loglik.1, loglik.2, df))
  
  if (loglik.1 < loglik.2) 
    return(1 - pchisq(2 * (loglik.2 - loglik.1), df))
  else return(1 - pchisq(2 * (loglik.1 - loglik.2), df))
  
}
