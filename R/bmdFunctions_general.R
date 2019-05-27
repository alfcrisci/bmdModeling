#' Main function to perform proast analysis
#' @param odt list as returned by f.scan()
#' @param shinyInput list with values assigned in the shiny app
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @param debug logical, if TRUE tryCatch errors and return current object ans.all; default value is FALSE
#' @return ans.all, list with all results that were obtained during the analysis
#' @export
f.proast <- function (odt = list(), shinyInput, track = FALSE, debug = FALSE) {
  
  if (track)	print("f.proast")
  
  if (debug) {
    
    tryCatch({
          
          ans.all <- f.ini(odt = odt, shinyInput = shinyInput, track = track)
          ans.all$user.nm <- deparse(substitute(odt))
          
          ans.all <- f.execute(ans.all, track = track)
          
          if (ans.all$cont){
            
            ans.all <- f.con(ans.all, track = track)
            
          } else {
            
            ans.all <- f.cat(ans.all, track = track)
            
          }
          
          if (track) 
            cat("\n\n\n    f.proast  END \n\n\n")
          
          return(ans.all)
          
        }, error = function(err){
          
          print(err)
          return(ans.all)
          
        })
    
  } else {
    
    ans.all <- f.ini(odt = odt, shinyInput = shinyInput, track = track)
    ans.all$user.nm <- deparse(substitute(odt))
    
    ans.all <- f.execute(ans.all, track = track)
    
    if (ans.all$cont){
      
      ans.all <- f.con(ans.all, track = track)
      
    } else {
      
      ans.all <- f.cat(ans.all, track = track)
      
    }
    
    if (track) 
      cat("\n\n\n    f.proast  END \n\n\n")
    
    return(ans.all)
    
    
  }
  
}

#' Define default parameter values needed in the proast functions
#' @param odt list as returned by f.scan()
#' @param shinyInput list with values assigned in the shiny app
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return updated version of ans.all, list with all results that were obtained
#' during the analysis
#' @export
f.ini <- function (odt = NULL, shinyInput, track = FALSE) {
  
  if (track)	print("f.ini")
  
  ans.all <- shinyInput
  
  ans.all$fit.ans <- 1  # type of error in response
  
  ans.all$group <- 0
  ans.all$model.type <- ifelse(ans.all$dtype == 3, 2, 1)
  ans.all$eps <- 1e-12
  ans.all$response.matr <- 0
  if(is.null(ans.all$ces.ans))
    ans.all$ces.ans <- 3
  
  
  # Boolean, set to TRUE if the likelihood is missing and start values should be adapted
  ans.all$errorAdjustStartValues <- FALSE
  
  
  # for quantal model
  if(!ans.all$cont){
    ans.all$kk <- NA	
    ans.all$nn.tot <- 0
    ans.all$kk.tot <- 0
    ans.all$nth <- 1
    ans.all$pi.full <- NA
    ans.all$th.par <- NA
    ans.all$sig.par <- NA
    ans.all$gr.txt <- character(0)
    ans.all$alfa.start <- 10
    ans.all$covariate <- 1
    ans.all$constr <- -Inf
    ans.all$shift <- 0
    ans.all$decr.zz <- TRUE
    ans.all$down <- FALSE
    ans.all$sp.fact <- 1
    ans.all$cex.1 <- NA
    ans.all$CI.plt <- FALSE
    ans.all$fct1.full <- NA
    ans.all$regr.start <- numeric()
    ans.all$ans.scale <- 0
    ans.all$def.exc <- NA
    ans.all$alfa.length <- 0
    
    # from f.choose.model
    if(ans.all$dtype == 3){
      ans.all$CES <- 0
      ans.all$ces.ans <- 1
    }
    ans.all$CED.model <- ((ans.all$model.type == 1) & ans.all$model.ans %in% 15:26)|
        ((ans.all$model.type == 2) & (ans.all$model.ans %in% c(12:15, 22:25)))
    
#		ans.all$par.start <- NA
    
  }
  
  # required for binary model
  ans.all$conf.lev <- 0.9
  
  # Likelihood optimization parameters
  ans.all$control.lst <- f.control(3)
  
  # lower/upper bounds for contrained parameters (only for categorical model)
  if(is.null(ans.all$lb))
    ans.all$lb <- NA
  if(is.null(ans.all$ub))
    ans.all$ub <- NA
  
  if(is.null(ans.all$sf.x)){
    
    ans.all$sf.x <- 1
    
  } 
  
  if(ans.all$cont | (!ans.all$cont & ans.all$dtype == 3)){
    
    ans.all$model.name <- switch(as.character(ans.all$model.ans),
        
        "1" = "Null model: y = a",
        "11" = "Full model: y = group mean",
        "13" = "E3-CED: y = a*exp(bx^d)",
        "15" = "E5-CED: y = a * [c-(c-1)exp(-bx^d)]",
        "23" = "H3-CED: a * (1 - x^d/(b^d+x^d))",
        "25" = "H5-CED: a * (1 + (c-1)x^d/(b^d+x^d))"
    
    )
    
  } else {
    
    ans.all$model.name <- switch(as.character(ans.all$model.ans),
        
        "1" = "null model, y = a", 
        "14" = "full model: y = observed incidence",
        "16" = "two-stage in terms of BMD",
        "18" = "log-logistic in terms of BMD", 
        "19" = "Weibull in terms of BMD",
        "21" = "log-probit in terms of BMD",
        "24" = "gamma model in terms of BMD", 
        "25" = "probit model in terms of BMD",
        "26" = "logistic model in terms of BMD"
    ) 
    
  }
  
  
  # For the plots
  ans.all$xy.lim <- NA
  ans.all$color <- 1:100
  ans.all$heading <- ans.all$model.name
  
  
  # Save data information
  if (is.list(odt) && length(odt$nvar) == 1) {
    
    ans.all$varnames <- odt$varnames
    ans.all$nvar <- odt$nvar
    ans.all$odt <- odt$data
    ans.all$info <- odt$info
    
  }
  
  # Remove missing values & Check whether x and y are numeric vectors: in server.R
  ans.all$data.0 <- odt$data
  ans.all$data.0 <- ans.all$data.0[order(ans.all$data.0[, ans.all$xans]), ]
  
  
  allFactors <- c(ans.all$fct1.no, ans.all$fct2.no, ans.all$fct3.no, 
      ans.all$fct4.no, ans.all$fct5.no)
  allFactors <- allFactors[allFactors != 0]
  
  
  # from this part on, from execute function of proast 61.3
  # scale x
  ans.all$x <- as.numeric(ans.all$data.0[, ans.all$xans]) / ans.all$sf.x
  ans.all$x.leg <- ans.all$varnames[ans.all$xans]
  
  if (ans.all$sf.x != 1){
    
    ans.all$x.leg <- paste0(ans.all$x.leg, "/", ans.all$sf.x)
    
  } 
  
  # set y
  ans.all$y <- ans.all$data.0[, ans.all$yans]
  ans.all$y.leg <- ans.all$varnames[ans.all$yans]
  
  
  # Check for one level factors & provide warning
  if(length(allFactors) != 0){
    
    if (length(allFactors) == 1){
      
      oneLevel <- nlevels(ans.all$data.0[,allFactors]) == 1
      
    } else {
      
      oneLevel <- apply(ans.all$data.0[,allFactors], 2, 
          function(x) nlevels(x) == 1)
      
    }
    
    if(any(oneLevel)){
      
      parameters <- paste(c("a", "b", "var", "c", "d")[oneLevel], collapse = ",")
      warning("The factor you chose as covariate on parameter(s)", parameters, 
          "has only one level\n you might have selected a subgroup for this factor")
      
    }
    
  }
  
  ## Define factors for parameters
  # Partly copied from f.execute() and f.change.settings()
  ans.all$factor.name <- ""
  ans.all$fct1.txt <- ""
  ans.all$fct2.txt <- ""
  ans.all$fct3.txt <- ""
  ans.all$fct4.txt <- ""
  ans.all$fct5.txt <- ""
  
  ans.all$warningCovariates <- ""
  
  
  for(i in 1:5){
    
    if (is.null(ans.all[[paste0("fct", i, ".no")]]))
      ans.all[[paste0("fct", i, ".no")]] <- 0
    
  }
  
  
#  cat("\nATTENTION: parameter c does not occur in this model\n\n")
#  cat("\nATTENTION: parameter d does not occur in this model\n\n")
  if (ans.all$cont) {
    
    if (ans.all$model.ans %in% c(1, 13, 23)){
      
      if(any(ans.all$fct4.no > 0))
        ans.all$warningCovariates <- paste(ans.all$warningCovariates,
            "Parameter c does not occur in the chosen model.")
      ans.all$fct4.no <- 0
      
    }
    
    if (ans.all$model.ans == 1) {
      
      if(any(ans.all$fct2.no > 0))
        ans.all$warningCovariates <- paste(ans.all$warningCovariates,
            "Parameter b does not occur in the chosen model.")
      ans.all$fct2.no <- 0
      
      if(any(ans.all$fct5.no > 0))
        ans.all$warningCovariates <- paste(ans.all$warningCovariates,
            "Parameter d does not occur in the chosen model.")
      ans.all$fct5.no <- 0
      
    }
    
  } else {
    
    if (ans.all$model.ans %in% c(1, 18, 25, 26)){
      
      if (any(ans.all$fct4.no > 0))
        ans.all$warningCovariates <- paste(ans.all$warningCovariates,
            "Parameter c does not occur in the chosen model.")
      ans.all$fct4.no <- 0
      
    }
    
    if (any(ans.all$fct5.no > 0))
      ans.all$warningCovariates <- paste(ans.all$warningCovariates,
          "Parameter d does not occur in the chosen model.")
    ans.all$fct5.no <- 0
    
    if (ans.all$model.ans == 1) {
      
      if (any(ans.all$fct2.no > 0))
        ans.all$warningCovariates <- paste(ans.all$warningCovariates,
            "Parameter b does not occur in the chosen model.")
      ans.all$fct2.no <- 0    
      
    }
    
  }
  
  
  maxVar <- ncol(ans.all$data.0)
  
  for (i in 1:5) {
    
    # Detect and define interactions
    if (length(ans.all[[paste0("fct", i, ".no")]]) > 1) {
      
      newName <- paste0("interaction_", 
          paste(ans.all$varnames[ans.all[[paste0("fct", i, ".no")]]], collapse = "*"))
      
      if (newName %in% ans.all$varnames) {
        
        ans.all[[paste0("fct", i, ".no")]] <- which(newName == ans.all$varnames)[1]
        
      } else {
        
        ans.all$data.0[, (maxVar+1)] <- interaction(ans.all$data.0[, ans.all[[paste0("fct", i, ".no")]]])
        ans.all[[paste0("fct", i, ".no")]] <- maxVar + 1
        ans.all$varnames[maxVar + 1] <- newName
        maxVar <- maxVar + 1
        
      }
      
    }
    
    # Define fct#
    if (ans.all[[paste0("fct", i, ".no")]] == 0) {
      
      if(i == 3 & !ans.all$cont){
        
        ans.all[[paste0("fct", i)]] <- 1
        
      } else {
        
        ans.all[[paste0("fct", i)]] <- rep(1, nrow(ans.all$data.0))
        
      }
      
    } else {
      
      column <- ans.all[[paste0("fct", i, ".no")]]
      ans.all[[paste0("fct", i)]] <- as.numeric(factor(ans.all$data.0[, column]))
      ans.all[[paste0("fct", i, ".txt")]] <- levels(factor(ans.all$data.0[, column]))
      ans.all[[paste0("fct", i, ".name")]] <- ans.all$varnames[column]
      ans.all$factor.name <- paste0(ans.all$factor.name, " factor", i, ": ", ans.all$varnames[column])
      
      
      # Order of levelNames, see f.boot.con()
#      if (is.null(ans.all$levelNames)) {
#        
#        ans.all$levelNames <- paste0(ans.all[[paste0("fct", i, ".name")]], "_", 
#            ans.all[[paste0("fct", i, ".txt")]])
#        
#      } else {
#        
#        newNames <- paste0(ans.all[[paste0("fct", i, ".name")]], "_", 
#            ans.all[[paste0("fct", i, ".txt")]])
#        
#        if (any(newNames != ans.all$levelNames)) {
#         
#        tmpLevelNames <- expand.grid(newNames, ans.all$levelNames)
#        ans.all$levelNames <- apply(tmpLevelNames[ncol(tmpLevelNames):1], 1, 
#            function(x) paste(x, collapse = " & "))
#        
#      }
#                
#      }
      
      if (i == 2)   # Only factor for parameter b will result in several bmd values 
        ans.all$levelNames <- paste0(names(ans.all$data.0)[ans.all[[paste0("fct", i, ".no")]]], ".", 
            ans.all[[paste0("fct", i, ".txt")]])
      
    }
    
  }
  
  # when including covariates ans.plt is 0 else 1: different colors for groups
  if(any(c(ans.all$fct1.no, ans.all$fct2.no, ans.all$fct3.no, ans.all$fct4.no, ans.all$fct5.no) > 0))
    ans.all$ans.plt <- 0	else ans.all$ans.plt <- 1
  
  # included covariate for parameter d?
  if(ans.all$fct5.no > 0)
    ans.all$covar.dd <- 1	else ans.all$covar.dd <- 0
  
  
  if (track)	print("f.ini : END")
  
  return(ans.all)
  
}

#' Calculate values for the response variable given the response data type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.execute <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.execute")
  
  
  with(ans.all, {
        
        yy <- y        
        
        # test if covariate is the same for parameters (a, b)|(a, d)|(a, c)
        twice <- FALSE
        if (fct1.no != 0) 
          twice <- fct1.no == fct2.no
        if (!twice & fct1.no > 1) 
          twice <- fct1.no == fct5.no
        if (!twice & fct1.no > 1) 
          twice <- fct1.no == fct4.no
        
        
        if (cont) {
          
          # Automatic detection limit
          detlim <- 0
          if (min(y) == 0) 
            detlim <- min(y[y > 0])
          detlim.col <- nvar + 1
          data.0[, detlim.col] <- detlim
          Vdetlim <- data.0[, detlim.col]
          
          ans.all$detlim <- detlim
          ans.all$Vdetlim <- Vdetlim
          
          
          ans.all$low.y <- 0.98 * min(yy)
          ans.all$upp.y <- 1.02 * max(yy)
          
          nn <- 0
          
          
          if (dtype %in% c(10, 15, 250, 260)) {
            # Summary data
            
            sd <- data.0[, sans]
            nn <- data.0[, nans]
            
            if (sd.se == 2){
              
              sd <- sd * sqrt(nn)
              
            }               
            
            if (dtype %in% c(10, 15)) {
              # Log-scale
              
              y.mn <- y
              CV <- sd/y.mn
              mn.log <- log(y.mn/sqrt(1 + CV^2))
              yy <- exp(mn.log)
              ans.all$mn.log <- mn.log
              ans.all$sd2.log <- log(CV^2 + 1)
              
            } else if (dtype == 250) {
              # Orig scale
              
              y.mn <- y
              sd2 <- sd^2
              ans.all$mn.log <- y.mn
              ans.all$sd2.log <- sd2
              ans.all$yy <- y.mn
              
            }
            
          } else {
            
            sd <- NA
            
          }
        }
        
        if (dtype %in% c(4, 6, 84)) {
          # Quantal data
          
          nn <- data.0[, nans]
          if (any(y > nn)) {
            
            stop("Number of responses is larger than sample size")
            
          }
          
          yy <- y/nn
          
        }
        
        if (dtype %in% 2:3) {
          # Binary/ordinal
          
          nn <- 1
          y.original <- y
          scores.orig <- sort(unique(y))
          
          ymn <- tapply(y, x, mean)
          xmn <- tapply(x, x, mean)
          ymn <- as.numeric(ymn)
          
          if (length(xmn) > 1){
            if (var(xmn != 0) & var(ymn != 0)){ 
              if (cor(xmn, ymn) < 0) {
                
                warning("Proast assumes zero to represent normal, and higher scores abnormal
                        \nYour data seem to have the opposite direction and therefore will be reversed for analysis")
                
                scores.orig <- rev(scores.orig)
                
              }
            }
          } 
          
          score <- 0
          
          for (ii in 1:length(scores.orig)) {
            y[y.original == scores.orig[ii]] <- score
            score <- score + 1
          }
          
          scores.mtr <- cbind(scores.orig, levels(factor(y)))
          dimnames(scores.mtr) <- list(NULL, c("orig.scores", 
                  "temp.scores"))
          dum.ord <- sum(scores.mtr[, 1] == scores.mtr[, 2]) != 
              length(scores.mtr[, 1])
          if (dum.ord) {
            
            warning("The original scores have been transformed for analysis as follows", scores.mtr)
            
          }
        }
        
        # from f.adjust.saved
        ans.all$nr.aa <- max(fct1)
        ans.all$nr.bb <- max(fct2)
        if (cont) {
          
          ans.all$nr.var <- max(fct3)
          
        } else {
          
          ans.all$nr.var <- 0
          
        }
        ans.all$nr.cc <- max(fct4)
        ans.all$nr.dd <- max(fct5)
        
        if (!cont) {
          
          nr.gr.out <- f.nr.gr(ans.all$nr.aa, ans.all$nr.bb, twice)
          ans.all$nr.gr <- nr.gr.out[1]
          
        } else {
          
          ans.all$sd <- sd
          ans.all$x.mn <- mean(x)
          
        }
        
        ans.all$yy <- yy
        ans.all$nn <- nn
        
        if (dtype %in% c(4, 6, 84)) {
          
          ans.all$kk <- y
          ans.all$y <- yy
          # don't consider censor data for now
          ans.all$cens <- 0
          
        }
        
        ans.all$cens.up <- NA        
        ans.all$twice <- twice
        
        
        if (dtype == 3) {
          
          ans.all$scores.orig <- scores.orig
          ans.all$scores.mtr <- scores.mtr
          ans.all$nth <- max(y)
          
        }
        
        if (dtype == 4 && max(ans.all$fct3) > 1){
          
          nth <- max(as.numeric(fct3))
          
        } 
        
        if (dtype == 6){
          
          ans.all$alfa.length <- 1
          
        } else {
          
          ans.all$fct1.full <- NA
          ans.all$fct2.full <- NA
          
        }
        
        ans.all$dtype.0 <- dtype
        
        if (track) 
          print("f.execute:  END")
        
        return(ans.all)
        
        
      })
}







#' Minimization of the log-likelihood function 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.nlminb <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.nlminb")
  
  with(ans.all, {
        
        scale.dum <- abs(1/par.start)
        scale.dum[scale.dum == Inf] <- 1000
        if (dtype == 6)	scale.dum[1] <- 1
        
        count <- 0
        
        if (cont) {
          
          fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
              # lower/upper for parameters
              lower = lb, upper = ub, control = control.lst, 
              x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
              fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
              mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
              Vdetlim = Vdetlim, CES = CES, twice = twice, 
              cens.up = cens.up, 
              # lower/upper bound for theta (in objective function f.lik.con)
              lb = lb, ub = ub, 
              par.tmp = par.start, increase = increase, x.mn = x.mn)
          loglik.check <- f.lik.con(fit.res$par, x, y, dtype, 
              fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
              nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
              ttt = ttt, fct4 = fct4, fct5 = fct5, 
              cens.up = cens.up, par.tmp = NA, increase = increase, 
              x.mn = x.mn)
          
          if (is.na(fit.res$obj)) 
            fit.res$obj <- 1e-12
          if (is.na(loglik.check)) 
            fit.res$obj <- 1e-12
          else if (fit.res$obj != loglik.check) 
            fit.res$obj <- 1e-12
          
          while (is.na(fit.res$par[1]) || (fit.res$obj == 0) || 
              abs(fit.res$obj) == 1e-12 || !is.finite(fit.res$obj)) {
            
            warning("Model is refitted with other scaling parameter")
            scale.dum <- rep(0.5, length(par.start))
            scale.dum <- 2 * scale.dum
            count <- count + 1
            
            fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum,
                lower = lb, upper = ub, control = control.lst, 
                x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                Vdetlim = Vdetlim, CES = CES, twice = TRUE, cens.up = cens.up, 
                lb = lb, ub = ub, par.tmp = fit.res$par, increase = increase, 
                x.mn = x.mn)
            
            if (is.na(fit.res$obj)) 
              fit.res$obj <- 1e-12
            if (count > 50) {
              fit.res$obj <- -1e+10
              fit.res$par <- 0
            }
            
          }
          
          ans.all$MLE <- fit.res$par
          ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
          
        } else {
          
          fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
              # lower/upper bounds for parameters
              lower = lb, upper = ub, 
              control = control.lst, 
              x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
              fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
              nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
              # added for testing ordinal model
              model.type = model.type, decr.zz = decr.zz, fct3.ref = fct3.ref,
              fct3 = fct3, fct4 = fct4, fct5 = fct5, CES.cat = CES.cat,
              CES = CES, ttt = ttt, twice = twice, 
#              cens = cens, 
              x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
              alfa.length = alfa.length, ces.ans = ces.ans, 
              CES1 = CES1, CES2 = CES2, 
              nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot)
          
          while (is.na(fit.res$par[1]) | (fit.res$obj == 0) | 
              (fit.res$obj == 1e-12) | is.na(fit.res$obj)) {
            
            warning("Model ", model.ans, " is refitted with other scaling parameter")
            
            # TODO include in shiny UI?
            if (model.ans == 14 && model.type == 1 && dtype == 6) 
              par.start[1] <- ans.all$alpha.start 
#                  eval(parse(prompt = paste("give start value for alfa", "  > ")))
            
            scale.dum <- rep(0.5, length(par.start))
            scale.dum <- 2 * scale.dum
            count <- count + 1
            
            fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                lower = lb, upper = ub, control = control.lst, 
                x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                # added for testing ordinal model
                model.type = model.type, decr.zz = decr.zz, fct3.ref = fct3.ref,
                CES = CES, ttt = ttt, twice = twice, 
#                cens = cens, 
                fct3 = fct3, fct5 = fct5,  CES.cat = CES.cat,
                x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
                alfa.length = alfa.length, ces.ans = ces.ans, 
                CES1 = CES1, CES2 = CES2,
                nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot)
            
            if (count > 50) {
              fit.res$obj <- -1e+10
              fit.res$par <- 0
            }
          }
          
          ans.all$MLE <- fit.res$par
          
          if (model.type == 2) {
            par.lst <- f.split.par(ans.all$MLE, nrp, nth, dtype)
            ans.all$regr.par <- par.lst$regr.par
            ans.all$th.par <- par.lst$th.par
            ans.all$sig.par <- par.lst$sig.par
          }
          else ans.all$regr.par <- ans.all$MLE
          
        }
        
        ans.all$loglik <- round(-fit.res$objective, 2)
        
        if (!is.finite(ans.all$loglik)) 
          ans.all$loglik <- -1e+12
        
        message <- fit.res$message
        ans.all$converged <- f.converged(message, fit.res$conv, track = track)
        
        if (dtype == 6 && (model.ans == 14 & model.type == 1)) {
          
          dum <- list()
          dum$loglik <- ans.all$loglik
          dum$MLE <- ans.all$MLE
          ans.all$full.model <- dum
          alfa.start <- ans.all$MLE[1]
          nn.dum <- table(x, x)
          nn.dum <- nn.dum[nn.dum != 0]
          if (length(nn.dum) > 15 && alfa.start > 20) 
            alfa.start <- 5
          ans.all$alfa.start <- alfa.start
          
        }
        
        ans.all$fitted <- TRUE
        
        if (track) 
          print("f.nlminb : END")
        
        return(ans.all)
        
      })
}


#' Define control parameters for nlminb(), which minimizes the log-likelihood function 
#' @param level integer, one of \code{1:3} determines which control parameters
#' are used; default value is 3
#' @return list, with control parameters 
#' @export
f.control <- function (level = 3) {
  
  lst <- list()
  switch(level, {
        
        lst$eval.max <- 50
        lst$iter.max <- 40
        lst$rel.tol <- 0.001
        lst$x.tol <- 0.015
        lst$step.min <- 0.00022
        
      }, {
        
        lst$eval.max <- 100
        lst$iter.max <- 75
        lst$rel.tol <- 1e-06
        lst$x.tol <- 0.00015
        lst$step.min <- 2.2e-07
        
      }, {
        
        lst$eval.max <- 1000
        lst$iter.max <- 750
        
      })
  
  return(lst)
  
}




#' Determine convergence type of the minimization algorithm
#' @param message object as returned from nlminb()$message
#' @param conv.out object as returned from nlminb()$conv
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return boolean, whether the optimization ended with convergence or not;
#' with attribute "message" which is the input object message
#' @export
f.converged <- function (message, conv.out, track = FALSE) {
  
  if (track) 
    print("f.converged")
  
  converged <- FALSE
  
  if (is.null(mode(message))) {
    
    warning("Convergence message is NULL")
    
  } else {
    
    converged <- 1 - conv.out
    
  }
  
  if (track) 
    print("f.converged:  END")
  
  attr(converged, "message") <- message
  
  return(converged)
  
}

#' Check whether any of the MLEs for the parameters hit the lower/upper constraint 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.hit.constr <- function (ans.all, track = FALSE) {
  
  if (track) 	print("f.hit.constr")
  
  with(ans.all, {
        
        if (!(model.ans %in% c(11, 14) && model.type == 1)) {
          
#          toReport <- c(MLE == lb && MLE != ub, MLE == ub && MLE != lb)
#		 		            
#          if (any(toReport)) {
#            
#            warning("A parameter estimate was equal to the 
#                    lower/upper constraint:", MLE[toReport])
#            
#          }
          
          for(side in c("upper", "lower")){		
            
            lst <- switch(side, 
                'upper' = {MLE == lb && MLE != ub},
                'lower' = {MLE == ub && MLE != lb}
            )
            if (model.type == 2) 
              lst[length(lst)] <- FALSE
            if (model.type == 2 && dtype %in% c(4, 6)) 
              lst[length(lst) - 1] <- FALSE
            if (any(lst)) {
              warning("The following parameter estimate was equal to the", side, "constraint: ",
                  MLE[lst], ".")
            }
            
          }	
          
        }
        
        if (track)	print("f.hit.constr:   END")
        
#        return(NULL)
        
      })
}


#' Construct regression parameter matrix 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.pars <- function (ans.all, track = FALSE) {
  
  if (track) print("f.pars")
  
  with(ans.all, {
        
        nrp <- length(regr.par)
        regr.par.matr <- numeric()
        
        gr.txt <- character(0)
        
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        if (length(ans.all$nr.cc) == 0) 
          nr.cc <- 1
        if (length(ans.all$nr.dd) == 0) 
          nr.dd <- 1
        if (length(ans.all$nr.var) == 0) 
          nr.var <- 1
        nr.subgr <- max(nr.aa, nr.bb, nr.cc, nr.dd)
        
        if (nr.aa == 1) 
          fct1.txt <- rep("", nr.subgr)
        if (nr.bb == 1) 
          fct2.txt <- rep("", nr.subgr)
        if (nr.cc == 1) 
          fct4.txt <- rep("", nr.subgr)
        if (nr.dd == 1) 
          fct5.txt <- rep("", nr.subgr)
        if (identical(fct1, fct2)) 
          fct2.txt <- rep("", nr.subgr)
        if (identical(fct1, fct4)) 
          fct4.txt <- rep("", nr.subgr)
        if (identical(fct4, fct2)) 
          fct4.txt <- rep("", nr.subgr)
        
        if (all(c(nr.aa, nr.bb, nr.cc, nr.dd, nr.var) == 1)) {
          
          ans.all$regr.par.matr <- matrix(regr.par, nrow = 1)
          ans.all$nr.gr <- 1
          return(ans.all)
          
        }
        
        kk <- 1
        
        
        if (ces.ans %in% (1:3)) {
          if (model.type == 2 & (model.ans %in% c(12:15, 22:25, 46))) {
            if (nrp > (nr.aa + nr.bb)) 
              par.rest <- regr.par[(nr.aa + nr.bb + 1):nrp]
            else par.rest <- numeric()
            if (nr.aa > 1 & nr.bb == 1) 
              for (ii in 1:nr.aa) {
                par.tmp <- c(regr.par[1], regr.par[ii + 1], 
                    par.rest)
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
              }
            if (nr.aa == 1 & nr.bb > 1) 
              for (jj in 1:nr.bb) {
                par.tmp <- c(regr.par[1], regr.par[jj + 1], 
                    par.rest)
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
              }
            if (nr.aa > 1 & nr.bb > 1) 
              for (jj in 1:nr.bb) {
                par.tmp <- c(regr.par[jj], regr.par[jj + 
                            nr.aa], par.rest)
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
              }
            if (nr.aa == 1 & nr.bb == 1) 
              regr.par.matr <- matrix(regr.par, nrow = 1)
            ans.all$regr.par.matr <- regr.par.matr
            gr.txt <- character(0)
            for (jj in 1:nr.bb) {
              f1 <- fct1[fct2 == jj]
              f1.lev <- levels(factor(f1))
              for (ii.index in f1.lev) {
                ii <- as.numeric(ii.index)
                if (!twice) {
                  gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                      sep = "-")
                }
                if (twice) {
                  if (length(ans.all$covar.txt) > 0) 
                    fct1.txt <- covar.txt
                  gr.txt[kk] <- paste(fct1.txt[ii], sep = "-")
                }
                kk <- kk + 1
              }
            }
            if (track)	print("f.pars special: END ")
            return(ans.all)
          }
        }
        
        if (all(c(nr.aa, nr.bb, nr.cc, nr.dd) == 1)){
          
          regr.par.matr <- matrix(regr.par, nrow = 1)
          
        } else if (nr.dd == 1) {
          
          if ((dtype %in% c(1, 5)) && (model.ans == 11) && 
              (nr.var > 1)) 
            fct1 <- fct3
          
          par.tmp <- rep(NA, length(regr.par))
          dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
          gr.txt <- character(0)
          for (jj in 1:nr.bb) {
            f1 <- fct1[fct2 == jj]
            f1.lev <- levels(factor(f1))
            for (ii.index in f1.lev) {
              ii <- as.numeric(ii.index)
              f4 <- fct4[fct1 == ii & fct2 == jj]
              f4.lev <- levels(factor(f4))
              for (mm.index in f4.lev) {
                mm <- as.numeric(mm.index)
                # MV added this because of warning: In rbind(regr.par.matr, par.tmp) :
#  number of columns of result is not a multiple of vector length (arg 2)
                cValue <- regr.par[nr.aa + nr.bb + mm]
                if(is.na(cValue) & !cont) cValue <- regr.par[nr.aa + nr.bb + 1]
                par.tmp <- c(regr.par[ii], regr.par[nr.aa + jj], 
                    cValue, dd)
                par.tmp <- par.tmp[!is.na(par.tmp)]
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
                gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                    fct4.txt[mm], sep = "-")
                kk <- kk + 1
              }
            }
          }
          
          if (!cont && model.ans == 14 && model.type == 1) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          } else if (cont && model.ans == 11) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          } else if(!cont){  #MV added this part for nr.cc = 2, but only 1 parameter for c
            n.col <- nrp - nr.aa - nr.bb - 1 + 3
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
          } else {
            n.col <- nrp - nr.aa - nr.bb - nr.cc + 3
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
          }
        }
        if (nr.cc == 1 && nr.dd > 1) {
          if ((dtype %in% c(1, 5)) && (model.ans == 11) && (nr.var > 1)) 
            fct1 <- fct3
          par.tmp <- rep(NA, length(regr.par))
          cc <- regr.par[nr.aa + nr.bb + 1]
          gr.txt <- character(0)
          for (jj in 1:nr.bb) {
            f1 <- fct1[fct2 == jj]
            f1.lev <- levels(factor(f1))
            for (ii.index in f1.lev) {
              ii <- as.numeric(ii.index)
              f5 <- fct5[fct1 == ii & fct2 == jj]
              f5.lev <- levels(factor(f5))
              for (mm.index in f5.lev) {
                mm <- as.numeric(mm.index)
                par.tmp <- c(regr.par[ii], regr.par[nr.aa + jj], cc, regr.par[nr.aa + nr.bb + 1 + mm])
                par.tmp <- par.tmp[!is.na(par.tmp)]
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
                gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], fct5.txt[mm], sep = "-")
                kk <- kk + 1
              }
            }
          }
          n.col <- nrp - nr.aa - nr.bb - nr.dd + 3
          if (model.ans == 14 && model.type == 1) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          } else regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
        }
        if (identical(fct1, fct2) && identical(fct2, fct4)) 
          if (nr.aa > 1 && nr.bb > 1 && nr.cc > 1) {
            ii <- 0
            jj <- nr.aa
            if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 41, 42)) 
              dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
            else dd <- NULL
            kk <- nr.aa + nr.bb
            zz <- 1
            gr.txt <- character(0)
            for (jj in 1:nr.cc) {
              if (nr.aa > 1) 
                ii <- ii + 1
              else ii <- 1
              if (nr.bb > 1) 
                jj <- jj + nr.aa
              else jj <- nr.aa + 1
              kk <- kk + 1
              # MV added this, see above
              cValue <- regr.par[kk]
              if(is.na(cValue) & !cont) cValue <- regr.par[nr.aa + nr.bb + 1]
              par.tmp <- c(regr.par[ii], regr.par[jj], cValue, dd)
              par.tmp <- par.tmp[!is.na(par.tmp)]
              regr.par.matr <- rbind(regr.par.matr, par.tmp)
              gr.txt[zz] <- paste(fct4.txt[zz])
              zz <- zz + 1
            }
          }
        if (identical(fct1, fct2) && identical(fct2, fct5)) {
          
          if (cont) {
#          # MV adapted this: proast code has regr.par.matr with too many rows
            if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 41, 42)) {
              regr.par.matr[,ncol(regr.par.matr)] <- regr.par[nr.aa + nr.bb + nr.cc + 1:nr.dd]
            } else {
              regr.par.matr[,ncol(regr.par.matr)] <- regr.par[nr.aa + nr.bb + 1:nr.dd]
            }
          } else {        
            
            if (nr.aa > 1 && nr.bb > 1 && nr.dd > 1) {
              ii <- 0
              jj <- nr.aa
              if (model.ans %in% c(3, 8, 13, 18, 23, 25, 41, 42)) {
                cc <- NULL
                kk <- nr.aa + nr.bb
              }
              if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 41, 42)) {
                cc <- regr.par[nr.aa + nr.bb + 1]
                kk <- nr.aa + nr.bb + 1
              }
              zz <- 1
              gr.txt <- character(0)
              for (jj in 1:nr.dd) {
                if (nr.aa > 1) 
                  ii <- ii + 1
                else ii <- 1
                if (nr.bb > 1) 
                  jj <- jj + nr.aa
                else jj <- nr.aa + 1
                kk <- kk + 1
                par.tmp <- c(regr.par[ii], regr.par[jj], cc, 
                    regr.par[kk])
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
                gr.txt[zz] <- paste(fct5.txt[zz])
                zz <- zz + 1
              }
            }
          }
        }
        ans.all$regr.par.matr <- regr.par.matr
        ans.all$gr.txt <- gr.txt
        ans.all$nr.gr <- nrow(regr.par.matr)
        
        if (track) 
          print("f.pars: END ")
        
        return(ans.all)
        
      })
}





#' Check whether the parameter value for c is too close to CES
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return boolean, TRUE if the the parameter value for c is too close to CES; 
#' if so warning message is printed 
#' @export
f.check.cc <- function (ans.all, track = FALSE) {
  
  with(ans.all, {
        
        if (track) 
          print("f.check.cc")
        
        cc.OK <- TRUE
        
        if (model.ans %in% c(15, 25)) {
          
          cc <- MLE[nr.var + nr.aa + nr.bb + 1]
          
          if (increase == 1 && (cc < 0.02 + (1 + CES))) {
            
            warning("The value of parameter c is too close to CES. 
                    This indicates (nonrandom) errors in the data or too high value of CES")
            cc.OK <- FALSE
            
          } else if (increase == -1 && (cc > -0.02 + (1 - abs(CES)))) {
            
            warning("The value of parameter c is too close to CES. 
                    This indicates (nonrandom) errors in the data or too high value of CES")
            cc.OK <- FALSE
            
          }
          
        }       
        
        if (track) 
          print("f.check.cc: END")        
        
        return(cc.OK)
        
      })
}


#' Wrapper function for calculating CI around parameter estimates, iterates 
#' f.profile.all()
#' @param ans.all list, with all results that were obtained during the analysis
#' @param nolog boolean; TRUE if calculations should be performed on the original
#' scale, FALSE if on the log-scale; default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.CI <- function (ans.all, nolog = FALSE, track = FALSE) {
  
  if (track)	print("f.CI")
  
  with(ans.all, {
        
        if (length(ans.all$MLE) > 25) {
          
          stop("Number of parameters exceeds 25")
          
        }
        
        profile.out <- f.profile.all(ans.all, nolog, track = track)
        
        updated <- FALSE
        count <- 0
        loglik.tmp <- profile.out$loglik
        
        while (sum(profile.out$MLE.new) != 0) {
          
          count <- count + 1
          
          if (track)	cat(" \n.... calculation of CI is re-started\n")
          
          ans.all$MLE <- profile.out$MLE.new
          
          if (cont)	ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
          
          ans.all$loglik <- profile.out$loglik
          
          if (ans.all$loglik < loglik.tmp - 0.02) {
            
            warning("Log-likelihood fluctuates, CI is not calculated")
            profile.out$conf.int <- matrix(NA, ncol = 2)
            profile.out$MLE.new <- 0
            
          } else {
            
            loglik.tmp <- ans.all$loglik
            profile.out <- f.profile.all(ans.all, nolog, track = track)
            updated <- TRUE
            
            if (count > 50) {
              
              warning("No global optimum found")
              profile.out$conf.int <- matrix(NA, ncol = 2)
              profile.out$MLE.new <- 0
              
            }
          }
        }
        
        if (updated & !cont) {
          
          switch(model.type, dum <- 1, dum <- 0)
          CED <- ans.all$MLE[(alfa.length + nr.aa + dum):
                  (alfa.length + nr.aa + dum + nr.bb)]
          CED.matr <- matrix(CED, ncol = 1)
          if (model.type == 1) 
            ans.all$regr.par <- ans.all$MLE[1:nrp]
          if (model.type == 2) {
            par.lst <- f.split.par(ans.all$MLE, nrp, nth, dtype)
            ans.all$regr.par <- c(1, par.lst$regr.par)
            ans.all$th.par <- par.lst$th.par
            ans.all$sig.par <- par.lst$sig.par
          }
          ans.all$CED <- CED
          ans.all$CED.matr <- CED.matr
          # MV removed according to proast62.10
#          ans.all$regr.par <- ans.all$MLE[1:nrp]
          
        }
        
        conf.int <- if (nolog)	profile.out$conf.int	else
              10^profile.out$conf.int
        
        # Replace missing values, if not both limits are missing
        conf.int <- t(apply(conf.int, 1, function(row){
                  isMissing <- is.na(row)
                  if(sum(isMissing) == 1){
                    row[isMissing] <- c(0, Inf)[isMissing]
                  }
                  return(row)
                }))
        
        ans.all$conf.int <- conf.int * sf.x
        ans.all$update <- updated
        ans.all$profile <- profile.out$profile
        
        if (track)	print("f.CI END")
        
        return(ans.all)
        
      })
}


#' Profile likelihood method to calculate confidence intervals 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param nolog boolean, whether the response is log-transformed or not;
#' default value is FALSE (ie log-transformed response)  
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.profile.all <- function (ans.all, nolog = FALSE, track = FALSE) {
  
  if (track) 
    print("f.profile.all")
  
  with(ans.all, {
        
        # Proast62.10 for probit model
        if (!cont & model.ans == 25) {
        
        max.runs <- 200
        count <- 0
        count.2 <- 0
        # Proast62.10
        CI.NA <- FALSE
        count.max <- 100
        crit <- 0.5 * qchisq(conf.lev, 1)
        jump.crit <- crit * 3
        small.step.low <- FALSE
        small.step.upp <- FALSE
        dist <- 1.01
        par.nr <- 0
        loglik.max <- loglik
        lb.orig <- lb
        ub.orig <- ub
        profile.out <- list()
        MLE.new <- 0
        tb <- "\t"
        
        sp.fact <- 1
        
        
        # Define number of groups
        if (cont) {
          
          sig2 <- mean(MLE[1:nr.var])
          
          if (group[1] == 0) {
            
            group <- nr.var + nr.aa + (1 : nr.bb)
            
          }
          
          if (nr.var > 1) 
            sp.fact <- fct3
          
          if (dtype == 5) {
            
            warning("Confidence intervals are calculated without accounting for \n
                    nested structure in the data (e.g. intralitter correlations) \n")
            
          }
          
        } else {
          
          if (dtype %in% c(4, 6)) 
            y <- kk/nn
          
          if (group[1] == 0) {
            
#            group <- nr.aa + (1 : nr.bb)
#            
#            if (dtype == 6) {
#              
#              group <- group + 1
#              
#            }
            from <- nr.aa + 1
            until <- nr.aa + nr.bb
            
            if (model.type == 2) {
              if (nr.aa == 1) {
                until <- until - nr.bb + 1
              }
              if (nr.aa > 1 & nr.bb == 1) {
                from <- 2
                until <- nr.aa + 1
              }
              if (nr.aa == 1 & nr.bb > 1) {
                from <- 2
                until <- nr.bb + 1
              }
              if (nr.aa > 1 & nr.bb > 1) {
                from <- nr.aa + 1
                until <- nr.aa + nr.bb
              }
            }
            
            if (dtype == 6) {
              from <- from + 1
              until <- until + 1
            }
            
            group <- from:until
            
          }
          
          sp.fact <- 1
          
        }
        
        if (nr.bb > 1) 
          sp.fact <- fct2
        if (nr.aa > 1) 
          sp.fact <- fct1
        
        conf.int <- matrix(NA, nrow = length(group), ncol = 2)
        
        for (jj in group) {
          
          par.nr <- par.nr + 1
          
          if (track) {
            
            if (model.type == 2 && nr.aa > 1 && nr.bb == 1) 
              cat("\ncalculating C.I.for group", fct1.txt[par.nr], 
                  "......\n")
            else if (nr.bb == 1 || quick.ans == 2) 
              cat("\nCalculating C.I.......\n")
            else cat("\nCalculating C.I.for group", fct2.txt[par.nr], 
                  "......\n")
            
          }
          
          # Proast62.10
          CED.upp.inf <- FALSE
          
          # Calculating lower limit
          if(track){
            cat("\n=========== lower limit =================\n")
            cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                        crit, 2), sep = tb)
          }
          
          lb <- lb.orig
          ub <- ub.orig
          lb[jj] <- MLE[jj]
          
          # Proast62.10
          if (lb[jj] < 0) 
            nolog <- T
          
          start <- dist * MLE
          loglik.low <- loglik.max
          CED.low <- MLE[jj]
          loglik.low.old <- numeric()
          CED.low.old <- numeric()
          # Proast62.10
          loglik.current <- loglik.max
          MLE.current <- MLE
          
          step.start <- if (cont && dtype %in% c(1, 5, 10))
                1.02 + 0.11 * sig2	else	1.08
          
          step <- step.start
          stop <- FALSE
          run <- 0
          
          while (!stop) {
            
            run <- run + 1
            lb[jj] <- lb[jj]/step
            ub[jj] <- lb[jj]
            start[jj] <- lb[jj]
            
            step.tmp <- step
            
            #Proast 61.3
#            if (lb[jj] < lb.orig[jj]) {
            
            # Proast 62.10
            while (lb[jj] < lb.orig[jj]) {
              step.tmp <- sqrt(step.tmp)
              lb[jj] <- lb.orig[jj] * step.tmp
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              step <- step.tmp
              
            }
            
            ans.all$par.start <- start
            ans.all$lb <- lb
            ans.all$ub <- ub
            
            # Proast61.3
#            nlminb.out <- f.nlminb(ans.all, track = track)
#            MLE.current <- nlminb.out$MLE
#            loglik.current <- nlminb.out$loglik
#            
#            if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
#              
#              if(track)	cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
#              start <- dist * nlminb.out$MLE
#              # Proast 61.3
            ##              loglik.low <- -Inf
#              # Proast 62.10
#              loglik.current <- -Inf
#              
#            }
            
            # Proast62.10
            if (cont) 
              loglik.try <- -f.lik.con(ans.all$par.start, 
                  x, y, dtype, fct1, fct2, fct3, model.ans, 
                  mn.log, sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                  twice = twice, ttt = ttt, fct4 = fct4, 
                  fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                  increase = increase, x.mn = x.mn, track = track)
            else loglik.try <- 0
            if (!is.finite(loglik.try)) {
              lb[jj] <- lb[jj] * step
              step <- sqrt(step)
              cat("\ntemporary fitting problem\n")
              start <- start * runif(length(start), 0.9, 
                  1.1)
            }
            else {
              nlminb.out <- f.nlminb(ans.all, track = track)
              MLE.current <- nlminb.out$MLE
              loglik.current <- nlminb.out$loglik
              start <- MLE.current * dist
              if (!is.finite(loglik.current) || loglik.current <= 
                  -1e+10) {
                cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
                start <- dist * nlminb.out$MLE
                loglik.current <- -Inf
              }
              
              if (loglik.current > loglik.max + 0.03) {
                
                ans.all$par.start <- MLE.current * dist
                if(track)	cat("\nlocal optimum found....\n")
                
                profile.out <- list()
                profile.out$MLE.new <- nlminb.out$MLE
                profile.out$loglik <- loglik.current
                profile.out$conf.int <- conf.int
                return(profile.out)
                
              }
              
              if (length(loglik.low) <= 2 || run > max.runs) {
                
                loglik.low <- c(loglik.current, loglik.low)
                CED.low <- c(lb[jj], CED.low)
                
              } 
              
              if (length(loglik.low) > 2) {
                
                if (loglik.current > loglik.low[2]) {
                  
                  loglik.low <- c(loglik.current, loglik.low[-(1:2)])
                  CED.low <- c(lb[jj], CED.low[-(1:2)])
                  
                } else if (loglik.current > loglik.low[1]) {
                  
                  loglik.low <- c(loglik.current, loglik.low[-1])
                  CED.low <- c(lb[jj], CED.low[-1])
                  
                } else {
                  
                  loglik.low <- c(loglik.current, loglik.low)
                  CED.low <- c(lb[jj], CED.low)
                  
                }
                
              }
              
              if (loglik.current == -Inf) {
                
                loglik.low <- loglik.low[-1]
                CED.low <- CED.low[-1]
                
              }
              
              if (length(loglik.low) > 1) {
                
                if (abs(loglik.low[1] - loglik.low[2]) > 0.25 * crit) 
                  step <- sqrt(step)
                if (abs(loglik.low[1] - loglik.low[2]) < 0.1 * crit) 
                  step <- step^1.26
                
              }
              
              if (track) {
                cat(signif(MLE.current, 3), round(loglik.current, 
                        2), round(loglik.max - crit, 2), sep = tb)
                
              }
            }
            
            if (step < 1.00001) {
              
              stop <- TRUE
              small.step.low <- TRUE
              
            }
            
            if (CED.low[1] < 1e-20) 
              stop <- TRUE
            
            if (length(loglik.low) > 30) 
              if (max(abs(diff(loglik.low[1:10]))) < 0.001) {
                stop <- TRUE
              }
            
            if (loglik.max - loglik.current > crit) {
              if (length(unique(loglik.low)) < 6) {
                if (track) 
                  cat("\n------- not enough points: ", length(loglik.low), 
                      " ------------------")
                
                CED.low <- c(CED.low, CED.low.old)
                loglik.low <- c(loglik.low, loglik.low.old)
                CED.low.old <- CED.low
                loglik.low.old <- loglik.low
                # Proast 61.3
#                power <- length(CED.low)/15
                # Proast 62.10
                power <- 1.1
                step <- step.start^power
                step.start <- step
                lb[jj] <- MLE[jj]
                start <- dist * MLE
                run <- 0
                count.2 <- count.2 + 1
                
                # Proast61.3
#                if (count.2 > 5) {
#                  
#                  count.2 <- -1
#                  stop <- TRUE
#                  loglik.low <- NA
#                  
#                }
                
                #Proast62.10
                if (count.2 > 10) {
                  cat("\n did not succeed to collect large enough number of points of the profile\n")
                  CI.NA <- TRUE
                  stop <- TRUE
                  loglik.low <- NA
                }
                
              } else if ((loglik.max - loglik.current) > jump.crit) {
                
                lb[jj] <- lb[jj] * step
                step <- step.start^0.5
                step.start <- step
                count <- count + 1
                
                if (step.start < 1.01) {
                  
                  CED.low[1] <- lb[jj]
                  loglik.low[1] <- loglik.current
                  count <- count.max
                  
                }
                
                if (count == count.max) 
                  stop <- TRUE
                run <- 0
                
              } else {
                
                loglik.diff <- diff(loglik.low[order(CED.low)])
                
#                print(loglik.low)
#                print(loglik.diff)
#                print(CED.low)
                
                if (sum(loglik.diff < -0.1) > length(CED.low) - 3) {
                  # missing values - error
                  
                  step.start <- 1 + (step.start - 1)/2
                  step <- step.start
                  run <- 0
                  
                } else {
                  
                  stop <- TRUE
                  
                }
              }
            }
            
            if (run > max.runs) {
              
              warning("Maximum number of runs reached in establishing profile")
              if (abs(loglik.low[1] - loglik.low[2]) < 0.01) 
                stop <- TRUE
            }
            
            jump.low <- FALSE
            if (count == count.max) {
              
              warning("The jump in the log-likelihood could not be mitigated")
              jump.low <- TRUE
              
            }
          }
          
          # Calculating upper limit
          if (track) {
            cat("\n\n=========== upper limit =================\n")
            cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                        crit, 2), sep = tb)
          }          
          
          if (MLE.new[1] != 0) 
            MLE <- MLE.new
          lb <- lb.orig
          ub <- ub.orig
          lb[jj] <- MLE[jj]
          start <- dist * MLE
          loglik.upp <- loglik.max
          CED.upp <- MLE[jj]
          CED.upp.old <- numeric()
          loglik.upp.old <- numeric()
          
          step.start <- if (cont && dtype %in% c(1, 5, 10))
                1.02 + 0.11 * sig2	else	1.08
          
          step <- step.start
          run <- 0
          stop <- FALSE
          
          count <- 0
          
          while (!stop) {
            
            run <- run + 1
            lb[jj] <- lb[jj] * step
            ub[jj] <- lb[jj]
            start[jj] <- lb[jj]
            
            
            step.tmp <- step
            #Proast 61.3
#            if (lb[jj] < lb.orig[jj]) {
#              
#              step.tmp <- sqrt(step.tmp)
#              lb[jj] <- lb.orig[jj] * step.tmp
#              lb[jj] <- lb[jj]/step.tmp
#              ub[jj] <- lb[jj]
#              start[jj] <- lb[jj]
#              step <- step.tmp
#              
#            }
            
            step.tmp <- step
            
            while (lb[jj] > ub.orig[jj]) {
              
              if(track)	cat("f.profile.all: value of parameter exceeds upper constraint, attempting to avoid this\n")
              step.tmp <- sqrt(step.tmp)
              lb[jj] <- lb[jj]/step.tmp
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              step <- step.tmp
              
            }
            
            ans.all$par.start <- start
            ans.all$lb <- lb
            ans.all$ub <- ub
            
            nlminb.out <- f.nlminb(ans.all, track = track)
            MLE.current <- nlminb.out$MLE
            loglik.current <- nlminb.out$loglik
            
            if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
              
              if(track)	cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
              start <- dist * nlminb.out$MLE
              loglik.upp <- -Inf
              
            }
            
            if (loglik.current > loglik.max + 0.03) {
              
              ans.all$par.start <- MLE.current * dist
              if(track)	cat("\nlocal optimum found....\n")
              
              profile.out <- list()
              profile.out$MLE.new <- nlminb.out$MLE
              profile.out$loglik <- loglik.current
              profile.out$conf.int <- conf.int
              return(profile.out)
              
            }
            
            if (length(loglik.upp) <= 2 || (run > max.runs)) {
              
              loglik.upp <- c(loglik.current, loglik.upp)
              CED.upp <- c(ub[jj], CED.upp)
              
            }
            
            if (length(loglik.upp) > 2) {
              
              if (loglik.current > loglik.upp[2]) {
                
                loglik.upp <- c(loglik.current, loglik.upp[-(1:2)])
                CED.upp <- c(ub[jj], CED.upp[-(1:2)])
                
              } else if (loglik.current > loglik.upp[1]) {
                
                loglik.upp <- c(loglik.current, loglik.upp[-1])
                CED.upp <- c(ub[jj], CED.upp[-1])
                
              } else if (loglik.current <= loglik.upp[1]) {
                
                loglik.upp <- c(loglik.current, loglik.upp)
                CED.upp <- c(ub[jj], CED.upp)
                
              }
            }
            
            if (loglik.current == -Inf) {
              
              loglik.upp <- loglik.upp[-1]
              CED.upp <- CED.upp[-1]
              
            }
            
            if (length(loglik.upp) > 1) {
              
              if (abs(loglik.upp[1] - loglik.upp[2]) > 0.25 * crit) 
                step <- sqrt(step)
              if (abs(loglik.upp[1] - loglik.upp[2]) < 0.1 * crit) 
                step <- step^1.26
            }
            
            if (track) {
              
              cat(signif(MLE.current, 3), round(loglik.current, 
                      2), round(loglik.max - crit, 2), sep = tb)
              
            }
            
            if (step < 1.00001) {
              
              stop <- TRUE
              small.step.upp <- TRUE
              
            }
            
            if (length(loglik.upp) > 30) 
              if (max(abs(diff(loglik.upp[1:10]))) < 0.001) 
                stop <- TRUE
            
            if (CED.upp[1] > 1e+20) 
              stop <- TRUE
            
            if (ub[jj] > 1e+50) {
              stop <- TRUE
            }
            
            if (loglik.max - loglik.current > crit) {
              if (length(unique(loglik.upp)) < 6) {
                if (track) 
                  cat("\n------- not enough points: ", length(loglik.upp), 
                      " ------------------")
                
                CED.upp <- c(CED.upp, CED.upp.old)
                loglik.upp <- c(loglik.upp, loglik.upp.old)
                CED.upp.old <- CED.upp
                loglik.upp.old <- loglik.upp
                # Proast 61.3
#                power <- length(CED.upp)/15
                # Proast 62.10
                power <- 0.9
                step <- step.start^power
                step.start <- step
                lb[jj] <- MLE[jj]
                start <- dist * MLE
                run <- 0
                count.2 <- count.2 + 1
                
                # Proast 61.3
#                if (count.2 > 5) {
#                  
#                  count.2 <- -1
#                  stop <- TRUE
#                  loglik.upp <- NA
#                  
#                }
                
                # Proast 62.10
                if (count.2 > 10) {
                  cat("\n did not succeed to collect large enough number of points of the profile\n")
                  CI.NA <- TRUE
                  stop <- TRUE
                  loglik.upp <- NA
                  CED.upp <- NA
                }
                
              } else if ((loglik.max - loglik.current) > jump.crit) {
                
                lb[jj] <- lb[jj]/step
                step <- step.start^0.5
                step.start <- step
                count <- count + 1
                
                if (step.start < 1.01) {
                  
                  CED.upp[1] <- lb[jj]
                  loglik.upp[1] <- loglik.current
                  count <- count.max
                  
                }
                
                if (count == count.max) 
                  stop <- TRUE
                
                run <- 0
                
              } else {
                
                loglik.diff <- diff(loglik.upp[order(CED.upp)])
                if (sum(loglik.diff < -0.1) > length(CED.upp) - 3) {
                  # missing values - error
                  
                  step.start <- 1 + (step.start - 1)/2
                  step <- step.start
                  run <- 0
                  
                }
                else {
                  
                  stop <- TRUE
                  
                }
              }
            }
            
            if (run > max.runs) {
              
              warning("Maximum number of runs reached in establishing profile")
              if (abs(loglik.upp[1] - loglik.upp[2]) < 0.01) 
                stop <- TRUE
              
            }
            
          }
          
          jump.upp <- FALSE
          if (count == count.max) {
            
            warning("The jump in the log-likelihood could not be mitigated")
            jump.upp <- TRUE
            
          }
          
          # Added in proast62.10
          if (MLE[jj] < 0) {
            loglik.low.neg <- loglik.low
            CED.low.neg <- CED.low
            loglik.low <- loglik.upp
            CED.low <- CED.upp
            loglik.upp <- loglik.low.neg
            CED.upp <- CED.low.neg
          }
          
          # Modify lower bounds
          if (!nolog) 
            CED.low <- log10(CED.low)
          
          if (max(abs(diff(loglik.low))) > 0.01) {
            
            loglik.low <- loglik.low[order(CED.low)]
            CED.low <- sort(CED.low)
            spline.low <- spline(CED.low, loglik.low)
            
            if (min(spline.low$x) > min(CED.low)) {
              
              spline.low$x <- c(min(CED.low), spline.low$x)
              spline.low$y <- c(loglik.low[1], spline.low$y)
              
            }
            
            if (length(unique(spline.low$y)) == 1) {
              
              if(track)	cat("\n\nonly one point in lower spline\n")
              ci.low <- list(y = Inf)
              conf.int[par.nr, 1] <- ci.low$y
              
            } else {
              
              ci.low <- approx(spline.low$y, spline.low$x, 
                  xout = loglik.max - crit)
              conf.int[par.nr, 1] <- ci.low$y
              
            }
            
            if (loglik.low[1] > loglik.max - crit) 
              conf.int[par.nr, 1] <- -Inf
            
            if (small.step.low) 
              conf.int[par.nr, 1] <- (CED.low[1])
            
          } else {
            
            spline.low <- list()
            spline.low$x <- CED.low
            spline.low$y <- rep(loglik.low[1], length(CED.low))
            ci.low <- list()
            ci.low$y <- min(CED.low)
            
          }
          
          # Modify upper bounds
          if (!nolog) 
            CED.upp <- log10(CED.upp)
          
          if (max(abs(diff(loglik.upp))) > 0.01) {
            
            loglik.upp <- loglik.upp[order(CED.upp)]
            CED.upp <- sort(CED.upp)
            spline.upp <- spline(CED.upp, loglik.upp)
            
            if (max(spline.upp$x) < max(CED.upp)) {
              
              spline.upp$x <- c(spline.upp$x, max(CED.upp))
              spline.upp$y <- c(spline.upp$y, loglik.upp[length(loglik.vec)])
              
            }
            
            if (length(unique(spline.upp$y)) == 1) {
              
              if(track)	cat("\n\nonly one point in upper spline\n")
              ci.upp <- list(y = Inf)
              conf.int[par.nr, 2] <- ci.upp$y
              
            } else {
              
              ci.upp <- approx(spline.upp$y, spline.upp$x, 
                  xout = loglik.max - crit)
              conf.int[par.nr, 2] <- ci.upp$y
              
            }
            
            if (loglik.upp[length(loglik.upp)] > loglik.max - crit) 
              conf.int[par.nr, 2] <- Inf
            
            if (small.step.upp) 
              conf.int[par.nr, 2] <- log10(CED.upp[1])
            
            if (CED.upp.inf) 
              conf.int[par.nr, 2] <- Inf
            
          } else {
            
            spline.upp <- list()
            spline.upp$x <- CED.upp
            spline.upp$y <- rep(loglik.upp[1], length(CED.upp))
            ci.upp <- list()
            ci.upp$y <- max(CED.upp)
            
          }
          
          # Bind results
          CED.vec <- c(CED.low, CED.upp)
          loglik.vec <- c(loglik.low, loglik.upp)
          
          if (jump.low) 
            conf.int[par.nr, 1] <- CED.vec[1]
          if (jump.upp) 
            conf.int[par.nr, 2] <- CED.vec[length(CED.vec)]
          
        } # End for loop over groups
        
        if (max(abs(diff(loglik.low))) < 0.01) 
          conf.int[par.nr, 1] <- -Inf
        
        if (max(abs(diff(loglik.upp))) < 0.01) 
          conf.int[par.nr, 2] <- Inf
        
        #Proast61.3
#        if (count.2 == -1) {
#          
#          conf.int[par.nr, 1] <- NA
#          conf.int[par.nr, 2] <- NA
#          
#        }
        #Proast 62.10
        if (CI.NA) {
          cat("\nCI could not be established, log-likelihood profile did not change\n")
          conf.int[par.nr, 1] <- NA
          conf.int[par.nr, 2] <- NA
        }
        
        profile.out <- list(conf.int = conf.int, MLE.new = MLE.new, 
            loglik = nlminb.out$loglik, profile = cbind(CED.vec, loglik.vec))
        
        
        
        # Proast61.3 for all other models
        } else {
          
          max.runs <- 200
          count <- 0
          count.2 <- 0
          count.max <- 100
          crit <- 0.5 * qchisq(conf.lev, 1)
          jump.crit <- crit * 3
          small.step.low <- FALSE
          small.step.upp <- FALSE
          dist <- 1.01
          par.nr <- 0
          loglik.max <- loglik
          lb.orig <- lb
          ub.orig <- ub
          profile.out <- list()
          MLE.new <- 0
          tb <- "\t"
          
          sp.fact <- 1
          
          # Define number of groups
          if (cont) {
            
            sig2 <- mean(MLE[1:nr.var])
            
            if (group[1] == 0) {
              
              group <- nr.var + nr.aa + (1 : nr.bb)
              
            }
            
            if (nr.var > 1) 
              sp.fact <- fct3
            
            if (dtype == 5) {
              
              warning("Confidence intervals are calculated without accounting for \n
                      nested structure in the data (e.g. intralitter correlations) \n")
              
            }
            
          } else {
            
            if (dtype %in% c(4, 6)) 
              y <- kk/nn
            
            if (group[1] == 0) {
              
#            group <- nr.aa + (1 : nr.bb)
#            
#            if (dtype == 6) {
#              
#              group <- group + 1
#              
#            }
              from <- nr.aa + 1
              until <- nr.aa + nr.bb
              
              if (model.type == 2) {
                if (nr.aa == 1) {
                  until <- until - nr.bb + 1
                }
                if (nr.aa > 1 & nr.bb == 1) {
                  from <- 2
                  until <- nr.aa + 1
                }
                if (nr.aa == 1 & nr.bb > 1) {
                  from <- 2
                  until <- nr.bb + 1
                }
                if (nr.aa > 1 & nr.bb > 1) {
                  from <- nr.aa + 1
                  until <- nr.aa + nr.bb
                }
              }
              
              if (dtype == 6) {
                from <- from + 1
                until <- until + 1
              }
              
              group <- from:until
              
            }
            
            sp.fact <- 1
            
          }
          
          if (nr.bb > 1) 
            sp.fact <- fct2
          if (nr.aa > 1) 
            sp.fact <- fct1
          
          conf.int <- matrix(NA, nrow = length(group), ncol = 2)
          
          for (jj in group) {
            
            par.nr <- par.nr + 1
            
            if (track) {
              
              if (model.type == 2 && nr.aa > 1 && nr.bb == 1) 
                cat("\ncalculating C.I.for group", fct1.txt[par.nr], 
                    "......\n")
              else if (nr.bb == 1 || quick.ans == 2) 
                cat("\nCalculating C.I.......\n")
              else cat("\nCalculating C.I.for group", fct2.txt[par.nr], 
                    "......\n")
              
            }
            
            # Calculating lower limit
            if(track){
              cat("\n=========== lower limit =================\n")
              cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                          crit, 2), sep = tb)
            }
            
            CED.upp.inf <- FALSE
            lb <- lb.orig
            ub <- ub.orig
            lb[jj] <- MLE[jj]
            
            start <- dist * MLE
            loglik.low <- loglik.max
            CED.low <- MLE[jj]
            
            loglik.low.old <- numeric()
            CED.low.old <- numeric()
            
            
            step.start <- if (cont && dtype %in% c(1, 5, 10))
                  1.02 + 0.11 * sig2	else	1.08
            
            step <- step.start
            stop <- FALSE
            run <- 0
            
            while (!stop) {
              
              run <- run + 1
              lb[jj] <- lb[jj]/step
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              
              step.tmp <- step
              
              if (lb[jj] < lb.orig[jj]) {
                
                step.tmp <- sqrt(step.tmp)
                lb[jj] <- lb.orig[jj] * step.tmp
                ub[jj] <- lb[jj]
                start[jj] <- lb[jj]
                step <- step.tmp
                
              }
              
              ans.all$par.start <- start
              ans.all$lb <- lb
              ans.all$ub <- ub
              
              nlminb.out <- f.nlminb(ans.all, track = track)
              MLE.current <- nlminb.out$MLE
              loglik.current <- nlminb.out$loglik
              
              if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
                
                if(track)	cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
                start <- dist * nlminb.out$MLE
                loglik.low <- -Inf
                
              }
              
              if (loglik.current > loglik.max + 0.03) {
                
                ans.all$par.start <- MLE.current * dist
                if(track)	cat("\nlocal optimum found....\n")
                
                profile.out <- list()
                profile.out$MLE.new <- nlminb.out$MLE
                profile.out$loglik <- loglik.current
                profile.out$conf.int <- conf.int
                return(profile.out)
                
              }
              
              if (length(loglik.low) <= 2 || run > max.runs) {
                
                loglik.low <- c(loglik.current, loglik.low)
                CED.low <- c(lb[jj], CED.low)
                
              } 
              
              if (length(loglik.low) > 2) {
                
                if (loglik.current > loglik.low[2]) {
                  
                  loglik.low <- c(loglik.current, loglik.low[-(1:2)])
                  CED.low <- c(lb[jj], CED.low[-(1:2)])
                  
                } else if (loglik.current > loglik.low[1]) {
                  
                  loglik.low <- c(loglik.current, loglik.low[-1])
                  CED.low <- c(lb[jj], CED.low[-1])
                  
                } else {
                  
                  loglik.low <- c(loglik.current, loglik.low)
                  CED.low <- c(lb[jj], CED.low)
                  
                }
                
              }
              
              if (loglik.current == -Inf) {
                
                loglik.low <- loglik.low[-1]
                CED.low <- CED.low[-1]
                
              }
              
              if (length(loglik.low) > 1) {
                
                if (abs(loglik.low[1] - loglik.low[2]) > 0.25 * crit) 
                  step <- sqrt(step)
                if (abs(loglik.low[1] - loglik.low[2]) < 0.1 * crit) 
                  step <- step^1.26
                
              }
              
              if (track) {
                cat(signif(MLE.current, 3), round(loglik.current, 
                        2), round(loglik.max - crit, 2), sep = tb)
                
              }
              
              if (step < 1.00001) {
                
                stop <- TRUE
                small.step.low <- TRUE
                
              }
              
              if (CED.low[1] < 1e-20) 
                stop <- TRUE
              
              if (length(loglik.low) > 30) 
                if (max(abs(diff(loglik.low[1:10]))) < 0.001) {
                  stop <- TRUE
                }
              
              if (loglik.max - loglik.current > crit) {
                if (length(unique(loglik.low)) < 6) {
                  if (track) 
                    cat("\n------- not enough points: ", length(loglik.low), 
                        " ------------------")
                  
                  CED.low <- c(CED.low, CED.low.old)
                  loglik.low <- c(loglik.low, loglik.low.old)
                  CED.low.old <- CED.low
                  loglik.low.old <- loglik.low
                  power <- length(CED.low)/15
                  step <- step.start^power
                  step.start <- step
                  lb[jj] <- MLE[jj]
                  start <- dist * MLE
                  run <- 0
                  count.2 <- count.2 + 1
                  
                  if (count.2 > 5) {
                    
                    count.2 <- -1
                    stop <- TRUE
                    loglik.low <- NA
                    
                  }
                  
                } else if ((loglik.max - loglik.current) > jump.crit) {
                  
                  lb[jj] <- lb[jj] * step
                  step <- step.start^0.5
                  step.start <- step
                  count <- count + 1
                  
                  if (step.start < 1.01) {
                    
                    CED.low[1] <- lb[jj]
                    loglik.low[1] <- loglik.current
                    count <- count.max
                    
                  }
                  
                  if (count == count.max) 
                    stop <- TRUE
                  run <- 0
                  
                } else {
                  
                  loglik.diff <- diff(loglik.low[order(CED.low)])
                  
#                print(loglik.low)
#                print(loglik.diff)
#                print(CED.low)
                  
                  if (sum(loglik.diff < -0.1) > length(CED.low) - 3) {
                    # missing values - error
                    
                    step.start <- 1 + (step.start - 1)/2
                    step <- step.start
                    run <- 0
                    
                  } else {
                    
                    stop <- TRUE
                    
                  }
                }
              }
              
              if (run > max.runs) {
                
                warning("Maximum number of runs reached in establishing profile")
                if (abs(loglik.low[1] - loglik.low[2]) < 0.01) 
                  stop <- TRUE
              }
              
              jump.low <- FALSE
              if (count == count.max) {
                
                warning("The jump in the log-likelihood could not be mitigated")
                jump.low <- TRUE
                
              }
            }
            
            # Calculating upper limit
            if (track) {
              cat("\n\n=========== upper limit =================\n")
              cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                          crit, 2), sep = tb)
            }          
            
            if (MLE.new[1] != 0) 
              MLE <- MLE.new
            lb <- lb.orig
            ub <- ub.orig
            lb[jj] <- MLE[jj]
            start <- dist * MLE
            loglik.upp <- loglik.max
            CED.upp <- MLE[jj]
            CED.upp.old <- numeric()
            loglik.upp.old <- numeric()
            
            step.start <- if (cont && dtype %in% c(1, 5, 10))
                  1.02 + 0.11 * sig2	else	1.08
            
            step <- step.start
            run <- 0
            stop <- FALSE
            
            count <- 0
            
            while (!stop) {
              
              run <- run + 1
              lb[jj] <- lb[jj] * step
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              
              
              step.tmp <- step
              if (lb[jj] < lb.orig[jj]) {
                
                step.tmp <- sqrt(step.tmp)
                lb[jj] <- lb.orig[jj] * step.tmp
                ub[jj] <- lb[jj]
                start[jj] <- lb[jj]
                step <- step.tmp
                
              }
              
              step.tmp <- step
              
              while (lb[jj] > ub.orig[jj]) {
                
                if(track)	cat("f.profile.all: value of parameter exceeds upper constraint, attempting to avoid this\n")
                step.tmp <- sqrt(step.tmp)
                lb[jj] <- lb[jj]/step.tmp
                ub[jj] <- lb[jj]
                start[jj] <- lb[jj]
                step <- step.tmp
                
              }
              
              ans.all$par.start <- start
              ans.all$lb <- lb
              ans.all$ub <- ub
              
              nlminb.out <- f.nlminb(ans.all, track = track)
              MLE.current <- nlminb.out$MLE
              loglik.current <- nlminb.out$loglik
              
              if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
                
                if(track)	cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
                start <- dist * nlminb.out$MLE
                loglik.upp <- -Inf
                
              }
              
              if (loglik.current > loglik.max + 0.03) {
                
                ans.all$par.start <- MLE.current * dist
                if(track)	cat("\nlocal optimum found....\n")
                
                profile.out <- list()
                profile.out$MLE.new <- nlminb.out$MLE
                profile.out$loglik <- loglik.current
                profile.out$conf.int <- conf.int
                return(profile.out)
                
              }
              
              if (length(loglik.upp) <= 2 || (run > max.runs)) {
                
                loglik.upp <- c(loglik.current, loglik.upp)
                CED.upp <- c(ub[jj], CED.upp)
                
              }
              
              if (length(loglik.upp) > 2) {
                
                if (loglik.current > loglik.upp[2]) {
                  
                  loglik.upp <- c(loglik.current, loglik.upp[-(1:2)])
                  CED.upp <- c(ub[jj], CED.upp[-(1:2)])
                  
                } else if (loglik.current > loglik.upp[1]) {
                  
                  loglik.upp <- c(loglik.current, loglik.upp[-1])
                  CED.upp <- c(ub[jj], CED.upp[-1])
                  
                } else if (loglik.current <= loglik.upp[1]) {
                  
                  loglik.upp <- c(loglik.current, loglik.upp)
                  CED.upp <- c(ub[jj], CED.upp)
                  
                }
              }
              
              if (loglik.current == -Inf) {
                
                loglik.upp <- loglik.upp[-1]
                CED.upp <- CED.upp[-1]
                
              }
              
              if (length(loglik.upp) > 1) {
                
                if (abs(loglik.upp[1] - loglik.upp[2]) > 0.25 * crit) 
                  step <- sqrt(step)
                if (abs(loglik.upp[1] - loglik.upp[2]) < 0.1 * crit) 
                  step <- step^1.26
              }
              
              if (track) {
                
                cat(signif(MLE.current, 3), round(loglik.current, 
                        2), round(loglik.max - crit, 2), sep = tb)
                
              }
              
              if (step < 1.00001) {
                
                stop <- TRUE
                small.step.upp <- TRUE
                
              }
              
              if (length(loglik.upp) > 30) 
                if (max(abs(diff(loglik.upp[1:10]))) < 0.001) 
                  stop <- TRUE
              
              if (CED.upp[1] > 1e+20) 
                stop <- TRUE
              
              if (ub[jj] > 1e+50) {
                stop <- TRUE
              }
              
              if (loglik.max - loglik.current > crit) {
                if (length(unique(loglik.upp)) < 6) {
                  if (track) 
                    cat("\n------- not enough points: ", length(loglik.upp), 
                        " ------------------")
                  
                  CED.upp <- c(CED.upp, CED.upp.old)
                  loglik.upp <- c(loglik.upp, loglik.upp.old)
                  CED.upp.old <- CED.upp
                  loglik.upp.old <- loglik.upp
                  power <- length(CED.upp)/15
                  step <- step.start^power
                  step.start <- step
                  lb[jj] <- MLE[jj]
                  start <- dist * MLE
                  run <- 0
                  count.2 <- count.2 + 1
                  
                  if (count.2 > 5) {
                    
                    count.2 <- -1
                    stop <- TRUE
                    loglik.upp <- NA
                    
                  }
                  
                } else if ((loglik.max - loglik.current) > jump.crit) {
                  
                  lb[jj] <- lb[jj]/step
                  step <- step.start^0.5
                  step.start <- step
                  count <- count + 1
                  
                  if (step.start < 1.01) {
                    
                    CED.upp[1] <- lb[jj]
                    loglik.upp[1] <- loglik.current
                    count <- count.max
                    
                  }
                  
                  if (count == count.max) 
                    stop <- TRUE
                  
                  run <- 0
                  
                } else {
                  
                  loglik.diff <- diff(loglik.upp[order(CED.upp)])
                  if (sum(loglik.diff < -0.1) > length(CED.upp) - 3) {
                    # missing values - error
                    
                    step.start <- 1 + (step.start - 1)/2
                    step <- step.start
                    run <- 0
                    
                  }
                  else {
                    
                    stop <- TRUE
                    
                  }
                }
              }
              
              if (run > max.runs) {
                
                warning("Maximum number of runs reached in establishing profile")
                if (abs(loglik.upp[1] - loglik.upp[2]) < 0.01) 
                  stop <- TRUE
                
              }
              
            }
            
            jump.upp <- FALSE
            if (count == count.max) {
              
              warning("The jump in the log-likelihood could not be mitigated")
              jump.upp <- TRUE
              
            }
            
            
            # Modify lower bounds
            if (!nolog) 
              CED.low <- log10(CED.low)
            
            if (max(abs(diff(loglik.low))) > 0.01) {
              
              loglik.low <- loglik.low[order(CED.low)]
              CED.low <- sort(CED.low)
              spline.low <- spline(CED.low, loglik.low)
              
              if (min(spline.low$x) > min(CED.low)) {
                
                spline.low$x <- c(min(CED.low), spline.low$x)
                spline.low$y <- c(loglik.low[1], spline.low$y)
                
              }
              
              if (length(unique(spline.low$y)) == 1) {
                
                if(track)	cat("\n\nonly one point in lower spline\n")
                ci.low <- list(y = Inf)
                conf.int[par.nr, 1] <- ci.low$y
                
              } else {
                
                ci.low <- approx(spline.low$y, spline.low$x, 
                    xout = loglik.max - crit)
                conf.int[par.nr, 1] <- ci.low$y
                
              }
              
              if (loglik.low[1] > loglik.max - crit) 
                conf.int[par.nr, 1] <- -Inf
              
              if (small.step.low) 
                conf.int[par.nr, 1] <- (CED.low[1])
              
            } else {
              
              spline.low <- list()
              spline.low$x <- CED.low
              spline.low$y <- rep(loglik.low[1], length(CED.low))
              ci.low <- list()
              ci.low$y <- min(CED.low)
              
            }
            
            # Modify upper bounds
            if (!nolog) 
              CED.upp <- log10(CED.upp)
            
            if (max(abs(diff(loglik.upp))) > 0.01) {
              
              loglik.upp <- loglik.upp[order(CED.upp)]
              CED.upp <- sort(CED.upp)
              spline.upp <- spline(CED.upp, loglik.upp)
              
              if (max(spline.upp$x) < max(CED.upp)) {
                
                spline.upp$x <- c(spline.upp$x, max(CED.upp))
                spline.upp$y <- c(spline.upp$y, loglik.upp[length(loglik.vec)])
                
              }
              
              if (length(unique(spline.upp$y)) == 1) {
                
                if(track)	cat("\n\nonly one point in upper spline\n")
                ci.upp <- list(y = Inf)
                conf.int[par.nr, 2] <- ci.upp$y
                
              } else {
                
                ci.upp <- approx(spline.upp$y, spline.upp$x, 
                    xout = loglik.max - crit)
                conf.int[par.nr, 2] <- ci.upp$y
                
              }
              
              if (loglik.upp[length(loglik.upp)] > loglik.max - crit) 
                conf.int[par.nr, 2] <- Inf
              
              if (small.step.upp) 
                conf.int[par.nr, 2] <- log10(CED.upp[1])
              
              if (CED.upp.inf) 
                conf.int[par.nr, 2] <- Inf
              
            } else {
              
              spline.upp <- list()
              spline.upp$x <- CED.upp
              spline.upp$y <- rep(loglik.upp[1], length(CED.upp))
              ci.upp <- list()
              ci.upp$y <- max(CED.upp)
              
            }
            
            # Bind results
            CED.vec <- c(CED.low, CED.upp)
            loglik.vec <- c(loglik.low, loglik.upp)
            
            if (jump.low) 
              conf.int[par.nr, 1] <- CED.vec[1]
            if (jump.upp) 
              conf.int[par.nr, 2] <- CED.vec[length(CED.vec)]
            
            
          } # End for loop over groups
          
          
          if (max(abs(diff(loglik.low))) < 0.01) 
            conf.int[par.nr, 1] <- -Inf
          
          if (max(abs(diff(loglik.upp))) < 0.01) 
            conf.int[par.nr, 2] <- Inf
          
          if (count.2 == -1) {
            
            conf.int[par.nr, 1] <- NA
            conf.int[par.nr, 2] <- NA
            
          }
          
          profile.out <- list(conf.int = conf.int, MLE.new = MLE.new, 
              loglik = nlminb.out$loglik, profile = cbind(CED.vec, loglik.vec))
          
        }
        
        
        if (track)
          print("f.profile.all:  END")
        
        
        return(profile.out)
        
      })
}




#' Determine the number of group ?
#' @param nr.aa values?
#' @param nr.bb values?
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @return numeric with bumber of groups?
#' @export
f.nr.gr <- function (nr.aa, nr.bb, twice) 
{
#	ab.ans <- 1 # why used for?
  nr.gr <- 1
#	if (nr.aa == 1) 
#		ab.ans <- 1
  if (nr.aa > 1 && nr.bb == 1) {
    nr.gr <- nr.aa
  }
  if (nr.bb > 1 && nr.aa == 1) {
    nr.gr <- nr.bb
  }
  if (nr.aa > 1 && nr.bb > 1) {
    if (twice) 
      nr.gr <- nr.aa
    if (!twice) 
      nr.gr <- nr.aa * nr.bb
  }
  return(nr.gr)
}



#' Define character vector with all parameter names
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return character vector, contains all parameter names
#' @export
f.text.par <- function (ans.all, track = FALSE) 
{
  if (track) print("f.text.par")
  with(ans.all, {
        
        if (!cont) 
          nr.var <- 0
        if (nr.aa == 1) {
          fct1.txt <- NULL
          sep.a <- NULL
        }
        if (nr.bb == 1 && length(xans) == 1) {
          fct2.txt <- NULL
          sep.b <- NULL
        }
        if (nr.var == 1) {
          fct3.txt <- NULL
          sep.v <- NULL
        }
        if (nr.cc == 1) {
          fct4.txt <- NULL
          sep.c <- NULL
        }
        if (nr.dd == 1) {
          fct5.txt <- NULL
          sep.d <- NULL
        }
        if (!cont && model.type == 1) {
          switch(model.ans, 
              text.par <- paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  "d"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "dd"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), 
              text.par <- c(rep("aa", nr.aa), rep("bb", 
                      nr.bb)), 
              text.par <- c(rep("aa", nr.aa), rep("bb", 
                      nr.bb), "cc"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("bb", fct2.txt[1:nr.bb], 
                      sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              { #model.ans = 14: full model
                if (length(ans.all$xx.tot) == 0) xx.tot <- x
#                if (dtype %in% c(4, 6, 84)) 
                text.par <- paste("group", 1:length(xx.tot))
              }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                      sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                      sep = "-"), "c"), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-"), "c", "d"), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-"), "c"), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-"), "c"), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-"), "c", "dd"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("BMD", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                      sep = "-"), "c"), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                      fct2.txt[1:nr.bb], sep = "-"), "cc"), 
              { #model.ans = 25: probit model
                if (ces.ans == 1) text.par <- c(paste("ED50", 
                          fct1.txt[1:nr.aa], sep = "-"), paste("bb", 
                          fct2.txt[1:nr.bb], sep = "-"))
                if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                          fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                          fct2.txt[1:nr.bb], sep = "-"))
              }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "BMDratio"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"))
        }
        if (cont || model.type == 2) {
          switch(model.ans, text.par <- paste("a", fct1.txt[1:nr.aa], 
                  sep = "-"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                      fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                      fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- rep("GM", 
                  length(regr.par)), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                      fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                      fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("CED", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), "c", "BMDratio"), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                      fct4.txt[1:nr.cc], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-"), 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("bb1", fct2.txt[1:nr.bb], sep = "-"), 
                  "bb2", "cc1", "cc2", "dd"), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                      fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                      fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("CED", fct2.txt[1:nr.bb], 
                      sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                  paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                      fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), "c"), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("bb1", fct2.txt[1:nr.bb], 
                      sep = "-"), "bb2", "cc", "dd1", "dd2"), text.par <- c(paste("aa", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("bb1", 
                      fct2.txt[1:nr.bb], sep = "-"), "bb2", "cc1", 
                  "cc2", "dd1", "dd2"), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-"), "cc", "dd"), text.par <- c(paste("a", 
                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", "d"), , , text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                  paste("d", fct5.txt[1:nr.dd], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-")), 
              text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                  paste("b", fct2.txt[1:nr.bb], sep = "-"), "c1", 
                  "c2", "d"), {
                ced.txt <- paste("RPF", fct2.txt[1:nr.bb], 
                    sep = "-")
                ced.txt[ref.lev] <- paste("CED", fct2.txt[ref.lev], 
                    sep = "-")
                text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                        sep = "-"), ced.txt, "c", paste("d", fct5.txt[1:nr.dd], 
                        sep = "-"))
              }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                  "c", paste("d", fct5.txt[1:nr.dd], sep = "-")))
          if (model.ans == 6 && ans.m6.sd == 2) 
            text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                "b", "q", "d ")
          if (model.ans == 47) 
            text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                "GCED", "q", "d ")
          if (cont && fit.ans == 1) 
            text.par <- c(paste("var", fct3.txt[1:nr.var], 
                    sep = "-"), text.par)
          if (length(xans) > 1) 
            text.par <- c(text.par, "RPF")
          if (!cont) {
            if (model.ans %in% c(12:15, 22:25) && nr.aa > 
                1 && nr.bb == 1) 
              for (ii in 2:(nr.aa + 1)) text.par[ii] <- paste("CED", 
                    ii - 1, sep = "")
            if (dtype == 4) 
              text.th <- "th"
            else text.th <- paste("th-", 1:nth, sep = "")
            text.par <- c(text.par, text.th, "sigma")
            if (ces.ans == 4) {
              con.lst <- f.start.con(ans.all.tmp, adjust = F, 
                  fitted = F, quick = F, no.warning = F)
              text.par <- con.lst$text.par
              if (model.ans > 1 & model.ans != 11) 
                text.par <- text.par[2:length(text.par)]
            }
          }
        }
        if (dtype == 6) 
          text.par <- c("alfa", text.par)
        if (track)  print("f.text.par:  END")
        return(text.par)
      })
}


