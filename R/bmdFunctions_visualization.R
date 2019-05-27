#' Wrapper function for plotting the results of proast
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @param plot.type, integer, defines the type of plot; for cont choose 1 (y vs x), 
#' for dtype 2 or 3 choose one of 3 (y.up vs. x), 5 (cumulative y vs. x); 
#' default is 1 for cont and 5 else
#' @return list, updated version of ans.all
#' @export
f.plot.all <- function (ans.all, track = FALSE, plot.type = NA) {
  
  if (track)	print("f.plot.all")
  
  # MV added default plot type choice
  if(is.na(plot.type)){
    
    if(ans.all$cont){
      
      ans.all$plot.type <- 1
      
    } else {
      
      ans.all$plot.type <- 5
      
    }
    
  } else {
    
    ans.all$plot.type <- plot.type
    
  }
  
  
  if (is.na(ans.all$xy.lim[1]) || ans.all$xy.lim[1] == 0) {
    
    x <- ans.all$x
    dum.contr <- min(x[x > 0])/5
    
    ans.all$xy.lim <- c(dum.contr, min(x), max(x), min(ans.all$y), max(ans.all$y))
    
  }
  
  with(ans.all, {
        
        if (cont) {
          
          ans.all$y.lim <- f.plot.con(ans.all, track = track)
          
          if (model.ans != 11){
            
            f.lines.con(ans.all, track = track)
            
          } 
          
          if (!is.na(CED[1]) & model.ans != 16) {
            
            f.cedlines.con(ans.all, track = track)
            
          }
          
        } else {
          
          if (dtype %in% c(4, 6, 84) | (dtype == 2 & plot.type < 3)) {
            
            ans.all <- f.plot.frq(ans.all, track = track)
            
            if (model.ans != 14) 
              ans.all$gr.txt <- f.lines.frq(ans.all, track = track)
            
            f.cedlines.bin(ans.all, track = track)
            
          } else if(dtype %in% 2:3) {
            
            ## in f.mm4.cat, code without combi:
            ans.all$combi <- F
            ans.all$categ.ans <- 0
#			f.plot.all(ans.all)
            
#			ans.all.plt <- ans.all #f.model.bb(ans.all = ans.all)
            
            # Error when plotting lines. In Proast61.3, NOTE:  Full model may not be applicable for dtype = 2
            if (!(model.ans == 14 & dtype == 2)) 
              ans.all$y.lim <- f.lines.cat(ans.all, track = track)
            
            title(ans.all$model.name)#f.mtext(ans.all)
            
            # in f.cat: code with combi.ans:
            ##			ans.all$combi.ans <- 1
            ##			ans.all$categ.ans <- 0
            ##			ans.all$plot.type <- 5
            ##			f.plot.all(ans.all, sep = F)
            ##			f.graph.window(4) equivalent to:
#			par(mfcol = c(2, 2), mar = c(3.2, 3.7, 2, 1.5), cex.main = 1.2,
#				cex.sub = 0.9, cex.lab = 0.9, cex = 0.9)
#			ans.all.plt <- ans.all
#			ans.all.plt$modelname <- ""
#			ans.all.plt$plot.type <- 3
#			ans.all.plt$ced.lines <- TRUE
#			ans.all.plt$main <- ""
#			ans.all.plt$y.lim <- f.lines.cat(ans.all.plt, track = track)
#			title(model.name, col = color[1], cex.main = 0.9, font.main = 1)
            ##			f.mtext(ans.all.plt, track = track)
#			ans.all.plt$plot.type <- 9
#			ans.all.plt$modelname <- ""
#			ans.all.plt$ced.lines <- F
#			ans.all.plt$y.lim <- f.lines.cat(ans.all.plt, track = track)
#			
#			if (!is.na(CED.matr[1])) {
#				th.par.cum <- cumsum(th.par)
#				CED.vect <- CED.matr[1, ]
#				if (ans.scale == 1) {
#					CED.vect <- log10(CED.vect)
#					x <- log10(x[x > 0]/sf.x)
#				}
#				for (ii in 1:length(CED.vect)) 
#					lines(c(min(x), CED.vect[ii]), rep(th.par.cum[ii], 2), lty = 2)
#				y.leg.2 <- paste(y.leg, " (latent var.uu)")
#				if (dtype == 3 & nr.gr == 1) {
#					CED.tmp <- round(CED.vect, rep(4, length(CED.vect)))
#					CED.tmp <- CED.tmp[CED.tmp < max(x)]
#					CED.left <- c(CED.tmp[1]/100, CED.tmp)
#					CED.right <- c(CED.tmp, max(x))
#					CED.mid <- (CED.left + CED.right)/2
#					if (ans.scale == 1) 
#						CED.mid <- 10^(CED.mid)
#					uu <- f.expect.con(model.ans, CED.mid, regr.par, 
#						fct1 = 1, fct2 = 1, CES = CES, track = track)
#					zz <- logb(uu)
#					if (ans.scale == 1) 
#						CED.mid <- log10(CED.mid)
#					for (ii in 1:length(scores.orig))
#						text(CED.mid[ii], zz[ii], scores.orig[ii])
#				}
#			}
#			th.par.cum <- cumsum(th.par)
#			if (mode(CED.matr) == "NULL") 
#				CED.vect <- 0
#			else if (is.na(CED.matr[1])) 
#				CED.vect <- 0
#			else CED.vect <- CED.matr[1, ]
#			ans.all.plt$plot.type <- 5
            
          }
          
        }
        
        if (track) 
          print("f.plot.all : END")
        
        return(ans.all)
        
      })
}


#' Plot the predicted response values (points) and estimated CIs 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return vector, range of the CI bounds on the predicted response values 
#' @export
f.plot.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.plot.con")
  
  with(ans.all, {
        
        if (is.null(x.leg)) 
          x.leg <- ""
        if (is.null(y.leg)) 
          y.leg <- ""
        
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        nr.var <- max(fct3)
        
        twice <- FALSE
        
        if (sum(fct1 != fct2) == 0) 
          twice <- TRUE
        
        f.means <- function(xx, yy, dtype, sd2.log, plot.type, nn) {
          
          if (track) 
            print("f.means within f.plot.con")
          
          
          mean.y <- yy
          
          if (dtype %in% c(1, 5, 15, 101)){
            
            mean.y <- exp(tapply(log(yy), xx, mean))
            
          } else if (dtype == 26){
            
            mean.y <- tapply(sqrt(yy), xx, mean)^2
            
          } else if (dtype == 25){
            
            mean.y <- tapply(yy, xx, mean)
            
          } 
          
          if (plot.type %in% c(1, 2, 5, 10, 11)) 
            mean.y <- mean.y
          if (plot.type %in% c(3, 4, 6)) 
            mean.y <- log10(mean.y)
          if (plot.type %in% c(7, 8)) 
            mean.y <- sqrt(mean.y)
          
          
          f.var <- function(dtype, xx, yy, sd2.log, nn) {
            
            if (track) 
              print("f.var within f.plot.con")
            
            if (dtype %in% c(10, 15, 110, 250, 260)) {
              
              SS <- sum(sd2.log * (nn - 1))
              var.within <- SS/sum(nn - 1)
              df <- sum(nn - 1)
              
            } else {
              
              if (dtype %in% c(1, 5)) {
                
                yy <- log(yy)
                
              } else if (dtype == 26){
                
                yy <- sqrt(yy)
                
              } 
              
              mn <- tapply(yy, xx, mean)
              nn <- tapply(yy, xx, length)
              Vmn <- rep(mn[1], nn[1])
              
              
              if (dtype %in% c(1, 5)){
                
                for (ii in 2:length(mn)){
                  Vmn <- c(Vmn, rep(mn[ii], nn[ii]))
                }
              } 
              
              resid <- yy - Vmn
              df <- sum(nn - 1)
              var.within <- (sum(nn) - 1) * var(resid)/df
              
            }
            
            out.lst <- list(Vsem = sqrt(var.within/nn), df = df)
            return(out.lst)
            
          }
          
          var.out <- f.var(dtype, xx, yy, sd2.log, nn)
          Vsem <- var.out$Vsem
          df <- var.out$df
          Vconf <- qt(0.975, df) * Vsem
          
          if (dtype %in% c(1, 5, 10, 15)) {
            
            conf.L <- mean.y/exp(Vconf)
            conf.U <- mean.y * exp(Vconf)
            
          } else if (dtype %in% c(25, 250)) {
            
            conf.L <- mean.y - Vconf
            conf.U <- mean.y + Vconf
            
          } else {
            
            conf.L <- (sqrt(mean.y) - Vconf)^2
            conf.U <- (sqrt(mean.y) + Vconf)^2
            
          }
          
          y.lim.tmp <- range(c(conf.L, conf.U))
          
          if (plot.type %in% c(3, 4, 6)) {
            conf.L <- log10(conf.L)
            conf.U <- log10(conf.U)
          }
          if (plot.type %in% c(7, 8)) {
            conf.L <- sqrt(conf.L)
            conf.U <- sqrt(conf.U)
          }
          
          y.lim.tmp <- c(min(conf.L), max(conf.U))
          
          out.lst <- list(conf.L = conf.L, conf.U = conf.U, 
              y.lim.tmp = y.lim.tmp, mean.y = mean.y)
          
          return(out.lst)
          
        }
        
        
        x.plt <- x
        y.plt <- yy
        
        switch(as.character(plot.type), 
            '1' = { #1
              x.plt <- x
              y.plt <- yy              
            }, 
            '3' = {  #3
              x.plt <- x
              y.plt <- log10(yy)              
            }, 
            '5' = { #5
              y.plt <- yy
              x.plt <- sqrt(x)
              if (x.leg != "") xleg <- paste("sqrt-", x.leg, sep = "") else xleg <- ""              
            })
        
        
        if (dtype %in% c(10, 110, 250, 260)) {
          
          plt.mns <- 3
          
        } else {
          
          plt.mns <- 1
          
        }
        
        cex.1 <- 0.6
        cex.2 <- 1.5
        
        
        mark <- c(1:2, 4:25, 33:500)
        
#        if (nr.aa == 1 && nr.bb == 1 && nr.var == 1) {
#          
#          if (max(sp.fact) > 1) {
#            
#            fct2 <- sp.fact
#            nr.bb <- max(fct2)
#            cat("\n\nno covariates, but subgroups are plotted distinctly\n\n")
#            
#          }
#          
#        }
        
        calculationsPlot <- f.means(xx = x.plt, yy = y.plt, dtype = dtype, 
            sd2.log = sd2.log, plot.type = plot.type, nn = nn)
        
        yRange <- range(c(y.plt, calculationsPlot$y.lim.tmp))
        
        plot(x.plt, y.plt, main = heading, ylim = yRange, 
            xlab = x.leg, ylab = y.leg, type = "n", col = color[1])
        
        
        if (nr.aa == 1 & nr.bb == 1 & nr.var > 1) {
          
          nr.bb <- nr.var
          fct2 <- fct3
          
        }
        
        if (nr.aa > 1 & nr.bb > 1 && twice) {
          
          nr.aa <- 1
          fct1 <- rep(1, length(x))
          
        }
        
        shift.tmp <- 0
        shift.abs <- 0
        bbb <- (max(x.plt, na.rm = TRUE) - min(x.plt, na.rm = TRUE))/70
        
        zz <- 0
        for (jj in 1:nr.bb) {
          
          for (ii in 1:nr.aa) {
            
            x.part <- x.plt[fct1 == ii & fct2 == jj]
            y.part <- y.plt[fct1 == ii & fct2 == jj]
            y.tmp <- yy[fct1 == ii & fct2 == jj]
            
            if (length(y.part) > 0) {
              
              zz <- zz + 1
              if (zz > 500) 
                cat("\nATTENTION: number of subgroups too large for plotting\n\n")
              type.dum <- "p"
              lt.y <- 1
              
              if (model.ans == 11) {
                lt.y <- 2
                type.dum <- "b"
              }
              
              x.part <- x.part + shift.tmp
              
              if (plt.mns < 3) {
                
                points(x.part, y.part, col = color[zz], pch = mark[zz], 
                    cex = cex.1, type = "p", lty = lt.y)
                
              }
              
              mean.x <- NA
              
              # Plot mean 
              # Note: changed compared to original code, because 
              # Error in xy.coords(x, y) : 'x' and 'y' lengths differ
              if (dtype %in% c(1, 5, 15, 25, 26)) {
                
                mean.x <- as.numeric(tapply(x.part, x.part, mean))
                mean.y <- f.means(x.part, y.tmp, dtype, plot.type = plot.type)$mean.y
                
                points(mean.x, mean.y, col = color[zz], pch = mark[zz], 
                    cex = 1, type = type.dum, lty = lt.y)
                
              } else {
                
                mean.x <- x.part
                mean.y <- y.part
                
                points(mean.x, mean.y, col = color[zz], pch = mark[zz], 
                    cex = 1, type = type.dum, lty = lt.y)
                
              }
              
              
              if (dtype %in% c(10, 250, 260)) {
                
                sd2.log.part <- sd2.log[fct1 == ii & fct2 == jj]
                nn.part <- nn[fct1 == ii & fct2 == jj]
                
              }
              
              out.lst <- f.means(x.part, y.tmp, dtype, sd2.log = sd2.log.part, 
                  plot.type, nn = nn.part)
              
              conf.L.part <- out.lst$conf.L
              conf.U.part <- out.lst$conf.U
              
              # Plot confidence limits
              for (qq in 1:length(conf.L.part)) {
                
                mean.x[qq] <- mean.x[qq]
                lines(rep(mean.x[qq], 2), c(conf.L.part[qq], 
                        conf.U.part[qq]), col = color[zz])
                lines(c(mean.x[qq] - bbb, mean.x[qq] + 
                            bbb), rep(conf.L.part[qq], 2), col = color[zz])
                lines(c(mean.x[qq] - bbb, mean.x[qq] + 
                            bbb), rep(conf.U.part[qq], 2), col = color[zz])
                
              }
              
            }
            
            shift.tmp <- shift.tmp + shift.abs
            
          }
          
        }
        
        
        ans.all$y.lim <- range(y.plt)
        ans.all$y.lim.CI <- range(y.plt)
        
        if (track) 
          print("f.plot.con:  END")
        
        return(calculationsPlot$y.lim.tmp)
        
      })
}


#' Plot dotted lines for the determined CED value for continous model
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return NULL
#' @export
f.cedlines.con <- function (ans.all, track = FALSE) {
  
  if (track)	print("f.cedlines.con")
  
  with(ans.all, {
        
        nr.subgr <- nrow(regr.par.matr)
        
        dum.contr <- xy.lim[1]
        # MV changed this
#        low.y <- min(y.lim)
        low.y <- xy.lim[4]
        
        for (jj in 1:nr.subgr) {
          
          CED.0 <- unlist(CED[jj])
          CES.0 <- CES
          
          ES.0 <- f.expect.con(model.ans, CED.0, regr.par.matr[jj,], 
              fct1 = 1, fct2 = 1, fct5 = 1, CES = CES.0, 
              increase = increase, track = track)
          
          ES.y <- rep(ES.0, 2)
          ES.x <- c(min(x), CED.0)
          CED.y <- c(low.y, ES.0)
          CED.x <- rep(CED.0, 2)
          
          lines(ES.x, ES.y, lty = 2, lwd = 1)
          lines(CED.x, CED.y, lty = 2, lwd = 1)
          
        }
        
        if (track)	print("f.cedlines.con: END")
        
        return(NULL)
        
      })
}


#' Plot the curve of the estimated model 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return NULL
#' @export
f.lines.con <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.lines.con")
  
  with(ans.all, {
        
        f.lines.tmp <- function(ans.all.tmp) {
          
          with(ans.all.tmp, {
                
                nbins <- 1000
                nbins.0 <- nbins/10
                
                xline <- seq(from = min(x, na.rm = TRUE), 
                    to = max(x, na.rm = TRUE), length = nbins)
                
                twice <- FALSE
                
                expect <- f.expect.con(model.ans, xline, regr.par, 
                    fct1 = rep(1, length(xline)), fct2 = rep(1, 
                        length(xline)), fct3 = 1, CES = CES, twice = twice, 
                    increase = increase, x.mn = x.mn, track = track)
                
                lines(xline, expect, lty = 1, col = col.tmp)
                
                return(NULL)
                
              })
          
        }
        
        ans.all.tmp <- ans.all
        
        ans.all.tmp$dum.contr <- xy.lim[1]
        ans.all.tmp$min.x <- xy.lim[2]
        ans.all.tmp$max.x <- xy.lim[3]
        if ((nr.aa > 1 || nr.bb > 1) && (nr.cc > 1 || nr.dd > 1) && !twice) {
          
          warning("Combination of covariate on c or d and yet another covariate 
                  \ndoes not allow plotting of all individual curves")
          
        }
        
        if (ans.plt == 1) 
          regr.par.matr <- matrix(regr.par, nrow = 1)
        else regr.par.matr <- f.pars(ans.all.tmp)$regr.par.matr
        
        count <- 0
        
        if (nr.aa == 1 && nr.bb == 1 && nr.var > 1){
          
          for (ii in 1:nr.var) {
            
            count <- count + 1
            ans.all.tmp$regr.par <- regr.par.matr[1, ]
            ans.all.tmp$col.tmp <- color[count]
            f.lines.tmp(ans.all.tmp)
            
          }
          
        } else {
          
          kk.max <- nrow(regr.par.matr)
          par.bb.old <- NA
          par.bb.new <- 0
          
          for (kk in 1:kk.max) {
            
            ans.all.tmp$regr.par <- regr.par.matr[kk, ]
            if (model.ans > 1) 
              par.bb.new <- ans.all.tmp$regr.par[2]
            
            if (ans.plt == 1) {
              if (!identical(par.bb.new, par.bb.old) || (model.ans == 
                    1)) {
                count <- count + 1
                ans.all.tmp$col.tmp <- color[count]
              }
            }
            else {
              count <- count + 1
              ans.all.tmp$col.tmp <- color[count]
            }
            
            f.lines.tmp(ans.all.tmp)
            par.bb.old <- par.bb.new
            
          }
        }
        
        if (track) 
          print("f.lines.con : END")
        
        return(NULL)
        
      })
}

#' Plot the predicted response values (points) and estimated CIs 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return updated ans.all
#' @export
f.plot.frq <- function (ans.all, track = FALSE){
  
  if(track)	print("f.plot.frq")
  
  with(ans.all, {
        
        if (length(fct1) == 1) 
          fct1 <- rep(1, length(x))
        
        if (length(fct2) == 1) 
          fct2 <- rep(1, length(x))
        
        dum.contr <- xy.lim[1]
        y.lim <- c(xy.lim[4], xy.lim[5])
        x.lim <- c(xy.lim[2], xy.lim[3])
        
        if (down)	y <- 1 - y
        if (max(sp.fact) > 1)	fct2 <- sp.fact
        if (max(fct1) > 1 & max(fct2) > 1 & (sum(fct1 == fct2) == length(fct1))) 
          twice <- T
        shift.abs <- 0
        shift.tmp <- 0
        mark <- c(1:25, 32:255)
        
        x.lim.plt <- x.lim
        y.lim.plt <- y.lim
        if (shift > 0) shift.abs <- (max(x) - min(x))/shift
        
        fct1 <- as.numeric(fct1)
        fct2 <- as.numeric(fct2)
        dum <- max(fct1, fct2)
        
        plot(x, y, #ylim = y.lim, xlim = x.lim, #removed because error if NA
            main = heading, xlab = x.leg, ylab = y.leg, type = "n", 
            col = color[1])
        
        if (shift == 0) 
          bbb <- (max(x) - min(x))/70	else 
          bbb <- (max(x) - min(x))/shift
        
        if (max(as.numeric(fct3)) > 1)	fct1 <- fct3
        zz <- 0
        
        for (jj in 1:max(fct2)) {
          for (ii in 1:max(fct1)) {
            x.tmp <- x[fct1 == ii & fct2 == jj]
            if (!is.na(x.tmp[1])) {
              zz <- zz + 1
              x.tmp <- x.tmp + shift.tmp
              y.tmp <- y[fct1 == ii & fct2 == jj]
              n.tmp <- nn[fct1 == ii & fct2 == jj]
              if (is.na(cex.1)) {
                if (dtype == 6) 
                  cex.1 <- 0.6
                else cex.1 <- 1.2
              }
              points(x.tmp, y.tmp, col = color[zz], pch = mark[zz], cex = cex.1)
              
              if (dtype == 4 & model.ans == 14) 
                lines(x.tmp, y.tmp, col = color[zz], lty = 2)
              
              if (CI.plt) {
                k.tmp <- kk[fct1 == ii & fct2 == jj]
                k.tmp <- k.tmp[!is.na(k.tmp)]
                n.tmp <- n.tmp[!is.na(n.tmp)]
                L025 <- numeric()
                L975 <- numeric()
                for (qq in 1:length(k.tmp)) {
                  LL.tmp <- f.LL.bin(k.tmp[qq], n.tmp[qq])
                  L025[qq] <- LL.tmp[1]
                  L975[qq] <- LL.tmp[2]
                }
                if (dtype %in% c(4, 6, 84)) {
                  points(x.tmp, L025, pch = "-", cex = 1.2, col = color[zz])
                  points(x.tmp, L975, pch = "-", cex = 1.2, col = color[zz])
                  
                  for (qq in 1:length(L025)) {
                    lines(rep(x.tmp[qq], 2), c(L025[qq], L975[qq]), col = color[zz])
                    lines(c(x.tmp[qq] - bbb, x.tmp[qq] + bbb), rep(L025[qq], 2), col = color[qq])
                    lines(c(x.tmp[qq] - bbb, x.tmp[qq] + bbb), rep(L975[qq], 2), col = color[zz])
                  }
                }
                shift.tmp <- shift.tmp + shift.abs
              }
            }
          }
        }
        if (!is.na(fct1.full[1])) 
          if (dtype == 6) {
            if (is.null(ans.all$cex.2)) 
              cex.2 <- 1.5
            zz <- 0
            for (jj in 1:max(fct2.full)) {
              for (ii in 1:max(fct1.full)) {
                x.tmp <- x.full[fct1.full == ii & fct2.full == 
                        jj]
                if (!is.na(x.tmp[1])) {
                  zz <- zz + 1
                  x.tmp <- x.tmp + shift.tmp
                  y.tmp <- pi.full[fct1.full == ii & fct2.full == 
                          jj]
                  if (sum(is.na(pi.full)) == 0) {
                    if (model.ans == 14) 
                      type.dum <- "b"
                    else type.dum <- "p"
                    points(x.tmp, y.tmp, cex = cex.2, pch = mark[zz], 
                        col = color[zz], type = type.dum, lty = 2)
                  }
                }
              }
            }
          }
        
        if (track)	print("plot.frq: END")
        ans.all$y.lim <- y.lim
        ans.all$y.lim.plt <- y.lim.plt
        ans.all$shift.abs <- shift.abs
        
        return(ans.all)
        
      })
}

#' Plot the lines
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return ans.all$gr.txt
#' @export
f.lines.frq <- function (ans.all, track = FALSE){
  
  if (track)	print("f.lines.frq")
  
  with(ans.all, {
        
        shift.abs <- 0
        if (model.ans == 14)	return(invisible())
        
        x.high <- max(x)
        
        f.lines.tmp <- function(model.ans, x, regr.par, 
            th.par, sig.par, plot.type, dtype, l.ty = l.ty,  
            color, x.high, shift, CES = CES, 
            twice, pi.full, down) {
          
          if (track)	print("f.lines.tmp in f.lines.frq")
          sat.mod <- FALSE
          if (model.ans == 14)	sat.mod <- TRUE
          
          if (!is.na(pi.full[1]) & sat.mod) {
            
            if (dtype %in% c(4, 84)) 
              xline <- x
            f(dtype == 6)
            xline <- tapply(x, x, mean)
            expect <- pi.full
            expect <- expect[order(xline)]
            xline <- sort(xline)
            
          } else {
            
            xline <- seq(from = min(x), to = x.high, length = 1000)
            
            expect <- f.expect.bin(model.ans, x = xline, regr.par, 
                CES = CES, ttt = ttt, twice = twice, ces.ans = ces.ans, 
                CES1 = CES1, CES2 = CES2)
            expectdum <- f.expect.bin(model.ans, x = 0, regr.par, 
                CES = CES, ttt = ttt, twice = twice, ces.ans = ces.ans, 
                CES1 = CES1, CES2 = CES2)
            if (down) {
              expect <- 1 - expect
              expectdum <- 1 - expectdum
            }
            pipi <- f.expect.bin(model.ans, x = x, regr.par, 
                CES = CES, ttt = ttt, twice = twice, ces.ans = ces.ans, 
                CES1 = CES1, CES2 = CES2)
            
            if (sum(is.na(expect)) > 0)	return()
            
            xline <- xline + shift
            
          }
          
          lines(xline, expect, lty = l.ty, col = color)
          
        }
        
        if (max(as.numeric(fct3)) > 1) {
          ans.all <- f.pars.frq(ans.all, track = track)
          regr.par.matr <- ans.all$regr.par.matr
          th.par.vec <- ans.all$th.par.vec
          kk.max <- length(th.par.vec)
          
          # MV added this to prevent errors for model.ans = 1 with covariates
          if(kk.max != nrow(regr.par.matr))
            regr.par.matr <- t(regr.par.matr)
          
        } else {
          ans.all <- f.pars(ans.all, track = track)
          regr.par.matr <- ans.all$regr.par.matr
          # MV changed this to prevent showing too many lines
#          kk.max <- length(regr.par.matr[, 1])
          kk.max <- nrow(unique(regr.par.matr))
          th.par.vec <- rep(th.par, kk.max)
        }
        
        shift.tmp <- 0
        
        for (kk in 1:kk.max) {
          
          f.lines.tmp(model.ans, x, regr.par.matr[kk, ], 
              th.par = th.par.vec[kk], sig.par, plot.type, 
              dtype = dtype, l.ty = l.ty, 
              col = color[kk], x.high = x.high, 
              shift = shift.tmp, CES = CES, 
              twice = twice, pi.full = pi.full, 
              down = down)
          shift.tmp <- shift.tmp + shift.abs
        }
        if (track)	print("f.lines.frq : END")
        return(ans.all$gr.txt)
        
      })
  
}

#' ?
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return updated ans.all
#' @export
f.pars.frq <- function (ans.all, track = FALSE) 
{
  if (track)	print("f.pars.frq")
  
  with(ans.all, {
        
        if (max(fct1) == 1 & max(fct2) == 1 & max(as.numeric(fct3)) == 1) 
          regr.par.matr <- matrix(regr.par, nrow = 1)
        else {
          
          if (max(as.numeric(fct3)) == 1) 
            pars <- c(regr.par, th.par)
          else {
            th.lvm <- f.th.lvm(th.par, fct3, fct3.ref, track = track)
            pars <- c(regr.par, th.lvm)
            fct3 <- as.numeric(factor(fct3))
          }
          nrp <- length(regr.par)
          regr.par.matr <- numeric()
          nr.aa <- max(fct1)
          nr.bb <- max(fct2)
          nr.th <- max(as.numeric(fct3))
          if (!covar.dd) {
            if ((model.ans == 11) & (nr.th > 1)) 
              fct1 <- fct3
            par.tmp <- rep(NA, length(pars))
            if (nrp > (nr.aa + nr.bb)) 
              par.rest <- pars[(nr.aa + nr.bb + 1):nrp]
            else par.rest <- numeric()
            zz <- 1
            for (jj in 1:nr.bb) {
              f1 <- fct1[fct2 == jj]
              f1.lev <- levels(factor(f1))
              for (ii.index in f1.lev) {
                ii <- as.numeric(ii.index)
                if (max(fct3) > 1) {
                  f3 <- fct3[fct1 == ii & fct2 == jj]
                  f3.lev <- levels(factor(f3))
                  for (kk.index in f3.lev) {
                    kk <- as.numeric(kk.index)
                    if (!twice) 
                      gr.txt[zz] <- paste(fct2.txt[jj], fct3.txt[kk], sep = "-")
                    if (twice) 
                      gr.txt[zz] <- fct2.txt[jj]
                    par.tmp <- c(pars[ii], pars[nr.aa + jj], 
                        par.rest, pars[nr.aa + nr.bb + length(par.rest) + kk])
                    regr.par.matr <- rbind(regr.par.matr, par.tmp)
                    zz <- zz + 1
                  }
                }
                if (max(fct3) == 1) {
                  par.tmp <- c(pars[ii], pars[nr.aa + jj], par.rest, th.par)
                  regr.par.matr <- rbind(regr.par.matr, par.tmp)
                  if (!twice) 
                    gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                        sep = "-")
                  if (twice) 
                    gr.txt[kk] <- fct1.txt[ii]
                }
              }
            }
            n.col <- nrp - nr.aa - nr.bb + 2 + 1
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
          }
          if (covar.dd) {
            cc <- pars[nr.aa + nr.bb + 1]
            for (jj in 1:nr.bb) {
              par.tmp <- c(pars[jj], pars[nr.aa + jj], cc, 
                  pars[nr.aa + nr.bb + 1 + jj])
              regr.par.matr <- rbind(regr.par.matr, par.tmp)
            }
          }
        }
        npar <- length(regr.par.matr[1, ])
        th.par.vec <- regr.par.matr[, npar]
        regr.par.matr <- regr.par.matr[, 1:(npar - 1)]
        
        # MV added this to prevent errors
        if(!is.matrix(regr.par.matr))
          regr.par.matr <- matrix(regr.par.matr, nrow = 1)
        
        ans.all$regr.par.matr <- regr.par.matr
        ans.all$gr.txt <- gr.txt
        ans.all$th.par.vec <- th.par.vec
        
        if (track)	print("f.pars.frq:   END")
        return(ans.all)
        
      })
}


#' compute th.lvm?
#' @param th.par parameter theta
#' @param fct3 initial value for theta
#' @param fct3.ref theta reference?
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return th.lvm
#' @export
f.th.lvm <- function (th.par, fct3, fct3.ref, track = FALSE) {
  if (track)	print("f.th.lvm")
  nr.th <- length(levels(as.factor(fct3)))
  fct3.lev <- as.numeric(levels(as.factor(fct3)))
  th.lvm <- rep(th.par, nr.th)
  th.lvm <- th.lvm - (log(fct3.ref) - log(fct3.lev))
  if (track)	print("f.th.lvm: END")
  return(th.lvm)
}

#' Plot dotted lines for the determined CED value for categorical model
#' @param ans.all list, with all results that were obtained during the analysis
#' @param inputResponse numeric, value for the response at estimated CED;
#' default value is NA
#' @param inputCED numeric, value for the CED; default value is NA
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return NULL
#' @export
f.cedlines.bin <- function (ans.all, inputResponse = NA, inputCED = NA, 
    track = FALSE) {
  
  if (track)	print("f.cedlines.bin")
  
  with(ans.all, {
        
        # MV added these conditions to use plot for model averaging
        if(is.na(inputResponse))
          response <- response.matr else
          response <- inputResponse
        
        if(is.na(inputCED))
          CED <- CED.matr else
          CED <- inputCED
        
        dum.contr <- xy.lim[1]
        nr.lev <- max(max(fct1, fct2))
        ES.y <- matrix(response[1], ncol = 2)
        ES.x <- matrix(c(x[1], CED[1]), ncol = 2)
        
        if (nr.lev > 1) 
          for (ii in (2:nr.lev)) {
            ES.y <- rbind(ES.y, rep(response[ii], 2))
            ES.x <- rbind(ES.x, c(x[1], CED[ii]))
          }
        
        CED.y <- matrix(c(min(y), response[1]), ncol = 2)
        CED.x <- matrix(rep(CED[1], 2), ncol = 2)
        
        if (nr.lev > 1) 
          for (ii in (2:nr.lev)) {
            CED.y <- rbind(CED.y, c(min(y), response[ii]))
            CED.x <- rbind(CED.x, rep(CED[ii], 2))
            
          }
        
        for (ii in (1:nr.lev)) {
          if (down) {
            ES.y[ii, ] <- 1 - ES.y[ii, ]
            CED.y[ii, 2] <- 1 - CED.y[ii, 2]
          }
          lines(ES.x[ii, ], ES.y[ii, ], lty = 2)
          lines(CED.x[ii, ], CED.y[ii, ], lty = 2)
        }
        
        return(NULL)
        
      })
}


#' Calculate confidence intervals for each of the plotted points for quantal response
#' @param k integer, number of events in given group
#' @param n integer, number of observations in given group
#' @param n.pi ? 
#' @param conf ?
#' @param step ?
#' @return vector with lower and upper limit
#' @export
f.LL.bin <- function (k = 4, n = 20, n.pi = 100, conf = 0.9, step = 1.2) 
{
  compl <- 0
  alfa <- (1 - conf)/2
  if ((k/n) > 0.5) {
    compl <- T
    k <- n - k
  }
  pi.min <- (k/n)
  pi.max <- 1
  pi.range <- 0.999 * seq(pi.min, pi.max, (pi.max - pi.min)/n.pi)
  prob1 <- pbinom(k, n, pi.range)
  if (sum(pi.range < 0 | pi.range > 1) > 0) {
  }
  dum1 <- pi.range[prob1 <= alfa]
  while (length(dum1) > n.pi/2) {
    pi.min <- pi.min/step
    pi.max <- pi.max/step
    pi.range <- 0.999 * seq(pi.min, pi.max, (pi.max - pi.min)/n.pi)
    prob1 <- pbinom(k, n, pi.range)
    if (sum(pi.range < 0 | pi.range > 1) > 0) {
    }
    dum1 <- pi.range[prob1 <= alfa]
  }
  dum2 <- pi.range[prob1 >= alfa]
  P2 <- (max(prob1[prob1 <= alfa]))
  P1 <- (min(prob1[prob1 >= alfa]))
  L.up <- ((P1 - alfa) * (min(dum1) - max(dum2)))/(P1 - P2) + 
      max(dum2)
  pi.min <- (k/n)
  pi.max <- 1
  pi.range <- 0.999 * seq(pi.min, pi.max, (pi.max - pi.min)/n.pi)
  if (k == 0) 
    L.low <- 0
  else {
    prob2 <- 1 - pbinom(k - 1, n, pi.range)
    if (sum(pi.range < 0 | pi.range > 1) > 0) {
    }
    if (sum(pi.range < 0 | pi.range > 1) > 0) {
    }
    dum1 <- pi.range[prob2 <= alfa]
    while (length(dum1) < n.pi/2) {
      pi.max <- pi.max/step
      pi.min <- pi.min/step
      pi.range <- 0.999 * seq(pi.min, pi.max, (pi.max - 
                pi.min)/n.pi)
      prob2 <- 1 - pbinom((k - 1), n, pi.range)
      if (sum(pi.range < 0 | pi.range > 1) > 0) {
      }
      dum1 <- pi.range[prob2 <= alfa]
    }
    dum2 <- pi.range[prob2 >= alfa]
    P2 <- (max(prob2[prob2 <= alfa]))
    P1 <- (min(prob2[prob2 >= alfa]))
    L.low <- min(dum2) - ((P1 - alfa) * (min(dum2) - max(dum1)))/(P1 - 
          P2)
  }
  if (compl) {
    L.low.tmp <- L.low
    L.low <- 1 - L.up
    L.up <- 1 - L.low.tmp
  }
  return(c(L.low, L.up))
}

#' plot lines for categorical model
#' 
#' Note: used in the f.plot.all function, for dtype == 2, 3, 
#' and where plot.type is set to 3, 5, 9?
#' so only this part of the code is retained
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return c(0, y.max)
#' @export
f.lines.cat <- function (ans.all, track = FALSE){
  
  if (track)	print("f.lines.cat")
  if (length(ans.all$cex.tmp) == 0) cex.tmp <- 1
  
  with(ans.all, {
        
        dum.contr <- xy.lim[1]
        x.lim <- c(xy.lim[2], xy.lim[3])
        y.max <- xy.lim[5]
        
        
        f.ced.lines <- function(regr.par, nth, dum.contr, CED, 
            x, response, plot.type, upper.uu = upper.uu, col.fct) {
          
          if (track)	print("f.ced.lines in f.lines.cat")
          
          switch(as.character(plot.type), 
              
              '3' = for (k in 1:nth)
                lines(rep(CED[k], 2), c(0, response[k]), 
                    lty = 2, col = col.fct), 
              
              '9' = for (k in 1:nth)
                lines(rep(CED[k], 2), c(min(logb(regr.par[2])), th.par[k]), 
                    lty = 2, col = col.fct)
          
          )
          
          switch(as.character(plot.type), 
              
              # 3
              '3' = for (k in 1:nth) 
                lines(c(0, CED[k]), c(response[k], response[k]), 
                    lty = 2, col = col.fct), 
              
              # 9
              '9' = for (k in 1:nth) 
                lines(c(0, CED[k]), rep(th.par[k], 2), lty = 2, col = col.fct)
          
          )
          # nothing is done for plot.type == 5
          
        }
        
        f.lines.tmp <- function(x, y, regr.par, th.par, sig.par, 
            nth, dum.contr, plot.type, col.fct, categ.ans = categ.ans, 
            first.plot, model.ans, model.type, x.leg, y.leg, 
            x.high, upper.uu = upper.uu, CES.cat, CES, reverse = reverse, 
            def.exc, cex.tmp, show, dtype, y.max, x.lim, twice) {
          
          cat.txt <- c("category 1", "category 2", "category 3", 
              "category 4", "category 5", "category 6", "category 7", 
              "category 8")
          
          if (track)	print("f.lines.tmp in f.lines.cat")
          
          xline <- seq(from = dum.contr, to = x.high, length = 100)
          
          if (model.type == 1) {
            
            #MV error when plot for dtype = 2, so added first part using dtype = 2 instead of default dtype = 6
            if(dtype == 2){
              
              pipi <- matrix(f.expect.bin(model.ans, x, regr.par, 
                      CES = CES, ces.ans = ces.ans, dtype = dtype, track = track), ncol = 1)
              piline <- matrix(f.expect.bin(model.ans, xline, 
                      regr.par, CES = CES, ces.ans = ces.ans, dtype = dtype, track = track), ncol = 1)
              
            } else {
              pipi <- matrix(f.expect.bin(model.ans, x, regr.par, 
                      CES = CES, ces.ans = ces.ans, track = track), ncol = 1)
              piline <- matrix(f.expect.bin(model.ans, xline, 
                      regr.par, CES = CES, ces.ans = ces.ans, track = track), ncol = 1)
            }
            
          } else {
            pipi <- f.expect.cat(model.type, model.ans, x, 
                regr.par, th.par, sig.par, CES = CES, CES.cat = CES.cat, 
                ces.ans = ces.ans, dtype = dtype, twice = twice, track = track)
            piline <- f.expect.cat(model.type, model.ans, 
                xline, regr.par, th.par, sig.par, CES = CES, 
                CES.cat = CES.cat, ces.ans = ces.ans, dtype = dtype, 
                twice = twice, track = track)
          }
          
          switch(as.character(plot.type), 
              
              # 3	
              '3' = {
                
                pi.up <- matrix(0, length(xline), (nth + 1))
                if (first.plot == T)
                  plot(xline, pi.up[, 1], 
                      xlim = x.lim, ylim = c(0, 1), xlab = x.leg, 
                      ylab = paste(y.leg, "(fraction of scores)"), 
                      type = "n", col = col.fct)
                
                for (k in 1:nth) {
                  dum <- nth + 1 - k
                  pi.up[, dum] <- piline[, dum] + pi.up[, (dum + 1)]
                }
                
                for (k in 1:nth)
                  lines(xline, pi.up[, k], col = col.fct, lty = k)
                y.max <- 1
                x.max <- x.high
                first.plot <- F
              }, 
              
              '5' = {
                
                y.leg <- paste(y.leg, " (cumu)")
#					if (!is.na(def.exc[1])) {
#						de.lev <- levels(factor(def.exc))
#						reverse <- T
#						lst.def <- def.exc == de.lev[1]
#						x.def <- x[lst.def]
#						y.def <- y[lst.def]
#						yy.def <- matrix(, length(x.def), nth)
#						pipi.tmp <- numeric(sum(lst.def))
#						for (i in (1:nth)) {
#							tmp <- pipi[, i]
#							tmp <- tmp[lst.def]
#							pipi.tmp <- cbind.data.frame(pipi.tmp, tmp)
#						}
#						pipi.def <- pipi.tmp[, 2:(nth + 1)]
#						ycum.def <- matrix(, length(x.def), nth)
#						picum.def <- matrix(, length(x.def), nth)
#						for (k in 1:nth) {
#							yy.def[, k] <- 1 * (y.def == k)
#							if (reverse) {
#								yy.def[, k] <- rev(yy.def[, k])
#								pipi.def[, k] <- rev(pipi.def[, k])
#							}
#							ycum.def[, k] <- cumsum(yy.def[, k])
#							picum.def[, k] <- cumsum(pipi.def[, k])
#							if (reverse) {
#								ycum.def[, k] <- rev(ycum.def[, k])
#								picum.def[, k] <- rev(picum.def[, k])
#							}
#						}
#						lst.exc <- def.exc == de.lev[2]
#						x.exc <- x[lst.exc]
#						y.exc <- y[lst.exc]
#						yy.exc <- matrix(, length(x.exc), nth)
#						pipi.tmp <- numeric(sum(lst.exc))
#						for (i in (1:nth)) {
#							tmp <- pipi[, i]
#							tmp <- tmp[lst.exc]
#							pipi.tmp <- cbind.data.frame(pipi.tmp, tmp)
#						}
#						pipi.exc <- pipi.tmp[, 2:(nth + 1)]
#						ycum.exc <- matrix(, length(x.exc), nth)
#						picum.exc <- matrix(, length(x.exc), nth)
#						for (k in 1:nth) {
#							yy.exc[, k] <- 1 * (y.exc == k)
#							ycum.exc[, k] <- cumsum(yy.exc[, k])
#							picum.exc[, k] <- cumsum(pipi.exc[, k])
#						}
#						x <- c(x.def, x.exc)
#						ycum <- rbind(ycum.def, ycum.exc)
#						picum <- rbind(picum.def, picum.exc)
#					} else {
                yy <- matrix(, length(x), nth)
                ycum <- matrix(, length(x), nth)
                picum <- matrix(, length(x), nth)
                for (k in 1:nth) {
                  yy[, k] <- 1 * (y == k)
                  ycum[, k] <- cumsum(yy[, k])
                  picum[, k] <- cumsum(pipi[, k])
                }
#					}
                y.max <- max(ycum)
                if (categ.ans < 2) {
                  if (first.plot == T) {
                    y.lim <- c(0, y.max)
                    plot(x, ycum[, 1], xlim = x.lim, ylim = y.lim, 
                        xlab = x.leg, ylab = y.leg, col = 1, cex = cex.tmp, 
                        type = "n")
                  }
                  for (k in 1:nth) {
                    points(x, ycum[, k], pch = paste(k), cex = cex.tmp, 
                        col = col.fct)
                    lines(x, picum[, k], col = col.fct, lty = k)
                  }
                }
#					if (categ.ans == 2) {
#						f.graph.window(nth)
#						for (ii in 1:nth) {
#							plot(x, ycum[, ii], ylim = c(0, max(ycum[, 
#															ii], picum[, ii])), xlab = x.leg, ylab = y.leg, 
#									type = "n", col = 1, cex = 0.6)
#							points(x, ycum[, ii], pch = "*", cex = 1.5, 
#									col = 1)
#							lines(x, picum[, ii], col = 1, lty = k)
#							title(paste("\ncategory ", ii), cex = 0.6)
#						}
#						first.plot <- F
#					}
                x.max <- x.high
#					if (plot.type == 6) x.max <- 1.1 * (x.high)
                if (categ.ans < 2) first.plot <- F
                
              },
              
              # 9	
              '9' = {
                
                zz <- f.expect.cat(model.type = 1, model.ans, xline, 
                    regr.par, th.par, latent = T, ces.ans = ces.ans, 
                    dtype = dtype, twice = twice, track = track)
                
                y.leg <- paste(y.leg, " (log(latent var))")
                upper.uu <- max(upper.uu, exp(max(zz)))
                if (sum(is.na(zz)) > 0) return()
                if (first.plot == T)
                  plot(xline, zz, xlim = x.lim, 
                      ylim = c(min(zz), logb(upper.uu)), xlab = x.leg, 
                      ylab = y.leg, type = "n", col = 1)
                
                lines(xline, zz, col = col.fct, lty = 1)
                y.max <- max(zz)
                x.max <- x.high
                first.plot <- F
                
              })
          
          xy.lim.tmp <- c(dum.contr, NA, x.max, 0, y.max, first.plot)
          if (track)	print("f.lines.tmp in f.lines.cat: END")
          
          return(xy.lim.tmp)
          
        }
        
        nth <- length(th.par)
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        regr.par.matr <- f.pars(ans.all)$regr.par.matr
        first.plot <- T
        kk <- 0
        fct2.lev <- levels(factor(fct2))
        x.high <- max(x)
        upper.uu <- max(regr.par.matr[, 1])
        
        if (plot.type %in% 5:6) {
          
          # Laure added: sort data according to x 
          xOld <- x
          x <- x[order(xOld)]
          y <- y[order(xOld)]
          
          y.max <- 0
          for (jj in fct2.lev) {
            f1 <- fct1[fct2 == as.numeric(jj)]
            f1.lev <- as.numeric(levels(factor(f1)))
            for (ii in f1.lev) {
              y.tmp <- y[fct1 == ii & fct2 == jj]
              yy <- matrix(, length(y.tmp), nth)
              ycum <- matrix(, length(y.tmp), nth)
              for (k in 1:nth) {
                yy[, k] <- 1 * (y.tmp == k)
                ycum[, k] <- cumsum(yy[, k])
              }
              y.max <- max(y.max, max(ycum))
              y.max <- max(max(ycum))
            }
          }
        }
        if (plot.type %in% 5:6) {
          kk <- 0
          nr.gr <- length(regr.par.matr[, 1])
          for (jj in fct2.lev) {
            f1 <- fct1[fct2 == as.numeric(jj)]
            f1.lev <- levels(factor(f1))
            for (ii in f1.lev) {
              kk <- kk + 1
              if (plot.type != 5 & plot.type != 6) 
                x.leg <- paste(x.leg, "\ngroup", kk)
              xy.lim.tmp <- f.lines.tmp(
                  x[fct1 == ii & fct2 == jj], y[fct1 == ii & fct2 == jj], 
                  regr.par.matr[kk, ], 
                  th.par, sig.par, nth, dum.contr, plot.type, 
                  col.fct = color[kk], categ.ans = categ.ans, 
                  first.plot = first.plot, model.ans, model.type, 
                  x.leg, y.leg, x.high, upper.uu = upper.uu, 
                  CES.cat = CES.cat, CES = CES, reverse = reverse, 
                  def.exc = def.exc[fct1 == ii & fct2 == jj], 
                  cex.tmp = cex.tmp, show = show, dtype = dtype, 
                  y.max = y.max, x.lim = x.lim, twice = twice)
              first.plot <- xy.lim.tmp[6]
            }
            if ((nr.gr > 1) & (kk < nr.aa * nr.bb)) {
              max.xx <- xy.lim[3]
              max.yy <- xy.lim[5]
              if (plot.type == 6) 
                max.xx <- max.xx + (max.xx - log10(dum.contr))/200
              else max.xx <- max.xx + (max.xx - min(x))/20
            }
          }
        }else	for (jj in fct2.lev) {
            
            f1 <- fct1[fct2 == as.numeric(jj)]
            f1.lev <- levels(factor(f1))
            
            for (ii in f1.lev) {
              
              kk <- kk + 1
              xy.lim.tmp <- f.lines.tmp(
                  x = x[fct1 == ii & fct2 == jj], y = y[fct1 == ii & fct2 == jj], 
                  regr.par = regr.par.matr[kk, ], th.par, sig.par, nth, dum.contr, plot.type, 
                  col.fct = color[kk], categ.ans = 0, first.plot = first.plot, 
                  model.ans, model.type, x.leg, y.leg, x.high, 
                  upper.uu = upper.uu, CES.cat = CES.cat, CES = CES, 
                  reverse = reverse, def.exc = NA, #def.exc[fct1 == ii & fct2 == jj], 
                  cex.tmp = cex.tmp, show = show, 
                  dtype = dtype, y.max = y.max, x.lim = x.lim, 
                  twice = twice
              )
              
              y.max <- xy.lim.tmp[5]
              first.plot <- xy.lim.tmp[6]
              if (!is.na(CED.matr[1])) {
                f.ced.lines(regr.par, nth, dum.contr, CED.matr[kk, ], 
                    x, response = response.matr[kk, ], plot.type, 
                    upper.uu = upper.uu, col.fct = color[kk]
                )
              }
              
            }
            
          }
        
        if (track)	print("f.lines.cat: END")
        
        return(c(0, y.max))
        
      })
  
}