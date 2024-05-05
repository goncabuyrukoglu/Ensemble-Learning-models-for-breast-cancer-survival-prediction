
##------------------------------------------------------------
## compare RF-SRC to Cox regression
## Illustrates C-index and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("party", logical.return = TRUE)
    & library("randomForestSRC", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))
  
{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }
  
  ## data, formula specifications
  data(GBSG2)
  trainrows <-  sample(1:nrow(GBSG2), .7*nrow(GBSG2))
  
  # GBSG2 <- GBSG2 %>% mutate(event = factor(cens, 
  #                                         labels =  c("censored", "failed")))
  #table(GBSG2$event)
  #GBSG2 <- GBSG2[,-10]
  surv.f <- as.formula(Surv(time, cens) ~.)
  pec.f <- as.formula(Hist(time, cens) ~ 1)
  
  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = GBSG2, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, GBSG2, ntree = 1000)
  fitcforest <- pecCforest(surv.f, data = GBSG2, controls = cforest_classical(ntree = 1000))
  ### estimate conditional Kaplan-Meier curves
  bst <- cforest(Surv(time, cens) ~ ., data = GBSG2,
                 control = cforest_unbiased(ntree = 50))
  
  treeresponse(bst, newdata = GBSG2[1:2,], OOB = TRUE)
  ### if you can't resist to look at individual trees ...
  party:::prettytree(bst@ensemble[[1]], names(bst@data@get("input")))
  vi <- varimp(bst)
  ### compare variable importances and absolute z-statistics
  layout(matrix(1:2))
  barplot(vi)
  barplot(abs(summary(coxph(Surv(time, cens) ~ ., data = GBSG2))$coeff[,"z"])) ### looks more or less the same
 
  
  
  plot.variable(rfsrc.obj, surv.type = "rel.freq", partial = TRUE)
  vip(bst)#Variable importance for survival forests (cforest)
  
  rsf.cindex <- get.cindex(time = GBSG2$time, censoring = GBSG2$cens, predicted = rfsrc.obj$predicted.oob)
  
  calPlot(rfsrc.obj, data = GBSG2)
  
  newData <- data.frame(horTh=factor("no",levels=c("no", "yes")),
                        age=as.integer(c(24,53,73)),
                        menostat=factor(cbind("Pre", "Post", "Post"),levels=c("Post","Pre")),
                        tsize=as.integer(rep(25, 3)),
                        tgrade=factor("II",levels=c("I","II", "III")),
                        pnodes=as.integer(rep(3, 3)), 
                        progrec=as.integer(rep(32, 3)), estrec=as.integer(rep(36, 3)))
  print(newData)
  
  
  
  
  cox.obj2 <- coxph(surv.f, data = GBSG2[trainrows,], x = TRUE)
  rfsrc.obj2 <- rfsrc(surv.f, GBSG2[trainrows,], ntree = 1000)
  fitcforest2 <- pecCforest(surv.f, data = GBSG2[trainrows,], controls = cforest_classical(ntree = 1000))
  
  
  
  
  
  
  cforestImpPlot <- function(x) {
    cforest_importance <<- v <- varimp(x)
    dotchart(v[order(v)])
  }
  cforestImpPlot(bst)
  
  ########################################################################################
  #################predicted survival curves for newdata#############################
  ########################################################################################
  
  pcox <- predictSurvProb(cox.obj, newdata = newData, times = quantile(GBSG2$time))
  prsf <- predictSurvProb(rfsrc.obj, newdata = newData, times = quantile(GBSG2$time))
  pcf <- predictSurvProb(fitcforest, newdata = newData, times = quantile(GBSG2$time))
  psp <- round(cbind(pcox,prsf, pcf),2)
  newNames <- paste("newData", 1:3)
  mat <- cbind( newData$age, psp)
  colnames(mat) <- c( "Age", rep("Coxph", 5), rep("rsf",5), rep("cforest",5))
  mat1 <- data.frame(newNames, mat)
  
  print(mat1)
  
  
  par(mfrow = c(1, 3))
  lapply(1:3, function(x) {
        plotPredictSurvProb(cox.obj, newdata = newData[x, ], col = 1)
        plotPredictSurvProb(rfsrc.obj, newdata = newData[x, ], add = TRUE, col = 2)
        plotPredictSurvProb(fitcforest, newdata = newData[x, ], add = TRUE,col = 3)
        })
  legend(x = "bottom",inset = 0, legend = c("Cox PH",
                                   "randomSurvivalForest",
                                   "cforest"),
         lwd = 3, col=1:3,bty="n",cex=1.1)
  
  ########################################################################################
  #########################################################################################
  ########################################################################################
  
  
  
  
  ########################################################################################
  ##################pec curves for bootstrap and bootstrap .632+ estimates#############
  ########################################################################################
  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  ptm <- proc.time()
  prederror.GBSG2 <- pec(list(cox.obj,rfsrc.obj, fitcforest), data = GBSG2[-trainrows,], formula = pec.f,
                       splitMethod = "bootcv", B = 500, M=200)
  print(prederror.GBSG2, times = seq(8, 2659,50))
  plot(prederror.GBSG2)
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  
  
  set.seed(17743)
  ptm <- proc.time()
  prederror.GBSG2.632 <- pec(list(cox.obj,rfsrc.obj, fitcforest), data = GBSG2[-trainrows,], formula = pec.f,
                         splitMethod = "Boot632", B = 500, M=200)
  print(prederror.GBSG2.632, times = seq(8, 2659,50))
  plot(prederror.GBSG2.632)
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  
  par(mfrow=c(2,2))
  plot(prederror.GBSG2, legend=FALSE)
  legend("bottomright", legend=c("Reference", "Cox PH","RSF", "Cforest"),
         col=c("black", "red","green", "blue"), lty=1)
  title(main = "GBSG2 bootstrap estimates")
  plot(prederror.GBSG2.632, legend=FALSE)
  legend("bottomright", legend=c("Reference", "Cox PH","RSF", "Cforest"),
         col=c("black", "red","green", "blue"), lty=1)
  title(main = "GBSG2 bootstrap.632 estimates")
  ########################################################################################
  #########################################################################################
  ########################################################################################
  
  
  
  
  #########################################################################################
  ###########cindex##############################################################################
  ########################################################################################
  
  # compute the bootstrap-crossvalidation estimate of
  # the C-index at different time points #
  ptm <- proc.time()
  set.seed(142)
  bcvCindex <- pec::cindex(list("Cox"=cox.obj,
                                "RSF"=rfsrc.obj,
                                "Cforest"= fitcforest), formula=surv.f,
                           data=GBSG2,
                           splitMethod="Bootcv",
                           B = 500, M=200, eval.times=c(50,quantile(GBSG2$time)[-1]))
  bcvCindex.632 <- pec::cindex(list("Cox"=cox.obj,
                                    "RSF"=rfsrc.obj,
                                    "Cforest"= fitcforest), formula=surv.f,
                               data=GBSG2,
                               splitMethod="Boot632",
                               B = 500, M=200, eval.times=c(50,quantile(GBSG2$time)[-1]))
  print(bcvCindex)
  plot(bcvCindex, legend=FALSE)
  legend("topleft", legend=c("Cox PH","RSF", "Cforest"),
           col=c("black", "red","green"), lty=1)
  title(main = "GBSG2 bootstrap estimates")
    
  print(bcvCindex.632)
  plot(bcvCindex.632, legend=FALSE)
  legend("topleft", legend=c("Cox PH","RSF", "Cforest"),
         col=c("black", "red","green"), lty=1)
  title(main = "GBSG2 bootstrap .632 estimates")
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  ########################################################################################
  #########################################################################################
  ########################################################################################
  
  
  
  
  
  #Rkare <- R2(prederror.GBSG2,times=seq(8, 2659,50))

  ######calibration plot cizdirmeye calistim ama olmadi....
  #veriyi train diye bolup oyle model kurdum hata bulmaya calisiyorum
  set.seed(17743)
  prederror2 <- pec(list(cox.obj2,rfsrc.obj2, fitcforest2), data = GBSG2[trainrows,], formula = pec.f,
                         splitMethod = "bootcv", B = 50)
  print(prederror2)
  plot(prederror2)
  
  
  cf2=calPlot(list("Cox regression"=cox.obj2,"RSF"=rfsrc.obj2),
              time=3,
              type="survival",
              data=GBSG2[-trainrows,])
  print(cf2)
  plot(cf2)
  calPlot(cox.obj2,time=3,data=GBSG2[-trainrows],type="survival")
  
  
  # oo <- subsample(rfsrc.obj)
  # # take a delete-d-jackknife procedure for example
  # vimpCI <- extract.subsample(oo)$var.jk.sel.Z
  # # Confidence Intervals for VIMP
  # plot.subsample(oo)
  # # take the variable "Month" for example for partial plot
  # plot.variable(rfsrc.obj, xvar.names = "days", partial = TRUE)
   

  
  
  
  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, GBSG2)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(GBSG2), nrow(GBSG2), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, GBSG2[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(GBSG2$time[-train], GBSG2$cens[-train], predict(cox.obj, GBSG2[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}






## illustrates the various splitting rules
## illustrates event specific and non-event specific variable selection

  ## events are transplant (1) and death (2)
  GBSG2$id <- NULL
  
  ## modified Gray's weighted log-rank splitting
  ## (equivalent to cause=c(1,1) and splitrule="logrankCR")
  GBSG2.cr <- rfsrc(Surv(time, cens) ~ ., GBSG2)
  
  ## log-rank cause-1 specific splitting and targeted VIMP for cause 1
  GBSG2.log1 <- rfsrc(Surv(time, cens) ~ ., GBSG2, 
                    splitrule = "logrank",importance="permute")
  
  GBSG2.log <- rfsrc(Surv(time, cens) ~ ., GBSG2, 
                     importance="permute")
  
  ## log-rank cause-2 specific splitting and targeted VIMP for cause 2
  #GBSG2.log2 <- rfsrc(Surv(time, cens) ~ ., GBSG2, 
   #                 splitrule = "logrank", cause = c(0,1), importance = TRUE)
  
  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(GBSG2.cr)$order[, 1],
                         vimp = 100 * GBSG2.log1$importance)
  lap <- print(var.perf[order(var.perf$md), ])
  stargazer(cox.obj,summary=TRUE, type = "text", title = "var per", out="cox.obj.txt")
  
  
  plot(GBSG2.log1)
  plot(GBSG2.log)# bunu ekledim sonuclara
  

  
  
  ############print results
  A2 <- summary(cox.obj2)
  write.table(round(A2$coefficients,4), file = "cox.obj2.txt", sep = ",", quote = FALSE)
