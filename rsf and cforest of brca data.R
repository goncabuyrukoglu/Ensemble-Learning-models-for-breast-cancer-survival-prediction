library(readr)
brca_metabric_clinical_data <- read_tsv("/Users/goncamert/Desktop/Gonca research/comparison of coxph and random survival forest/metabric dataset and papers/brca_metabric_clinical_data.tsv")
View(brca_metabric_clinical_data)
br <- data.frame(brca_metabric_clinical_data)
library(dplyr)

###########data preprocessing##########
brnew <- br[complete.cases(br$Patient.s.Vital.Status),]
brnew$age <- as.integer(round(brnew$Age.at.Diagnosis))
brnew1 <- brnew %>% select(-Study.ID, -Sample.ID,-Integrative.Cluster, -Primary.Tumor.Laterality,-Oncotree.Code,  -Age.at.Diagnosis,-Type.of.Breast.Surgery, -Cancer.Type, -Pam50...Claudin.low.subtype, -Cancer.Type.Detailed, -Cohort, -ER.status.measured.by.IHC, -HER2.status.measured.by.SNP6, -Tumor.Other.Histologic.Subtype, -Number.of.Samples.Per.Patient,-Sample.Type, -Sex, -Patient.s.Vital.Status, -X3.Gene.classifier.subtype)
brnew1$Cellularity <-  factor(brnew1$Cellularity, levels = c("Low", "Moderate", "High"))
brnew1$Chemotherapy <- factor(brnew1$Chemotherapy, levels = c("NO", "YES"))
brnew1$ER.Status <- factor(brnew1$ER.Status, levels = c("Negative", "Positive"))
brnew1$HER2.Status <- factor(brnew1$HER2.Status, levels = c("Negative", "Positive"))
brnew1$Hormone.Therapy <- factor(brnew1$Hormone.Therapy, levels = c("NO", "YES"))
brnew1$Inferred.Menopausal.State <- factor(brnew1$Inferred.Menopausal.State, levels = c("Post", "Pre"))
brnew1$PR.Status <- factor(brnew1$PR.Status, levels = c("Negative", "Positive"))
brnew1$Radio.Therapy <- factor(brnew1$Radio.Therapy, levels = c("NO", "YES"))
brnew1$Tumor.Stage[brnew1$Tumor.Stage==4] <-3 
brnew1$Tumor.Stage[brnew1$Tumor.Stage==0] <-1
brnew1$Tumor.Stage <- factor(brnew1$Tumor.Stage)
brnew1$Relapse.Free.Status <- substr(brnew1$Relapse.Free.Status, 1,1)
brnew1$Relapse.Free.Status<- as.integer(brnew1$Relapse.Free.Status)
brnew1$Overall.Survival.Status <- substr(brnew1$Overall.Survival.Status, 1,1)
brnew1$Neoplasm.Histologic.Grade <- factor(brnew1$Neoplasm.Histologic.Grade)
#brnew1 <- subset(brnew1, select=-c(Overall.Survival.Status, Overall.Survival..Months.))
#####delete missing values 
brn <- na.omit(brnew1)
na_count <-sapply(brn, function(y) sum(length(which(is.na(y)))))

###pool the event 
ii <- brn$Relapse.Free.Status==1 & brn$Overall.Survival.Status==1 
io <- brn$Relapse.Free.Status==1 & brn$Overall.Survival.Status==0
oi <- brn$Relapse.Free.Status==0 & brn$Overall.Survival.Status==1
oo <- brn$Relapse.Free.Status==0 & brn$Overall.Survival.Status==0
brn$time <- rep(0, dim(brn)[1])
brn$time[ii] <- brn$Relapse.Free.Status..Months.[ii]
brn$time[io] <- brn$Relapse.Free.Status..Months.[io]
brn$time[oi] <- brn$Overall.Survival..Months.[oi]
brn$time[oo] <- brn$Overall.Survival..Months.[oo]
brn$time <- 30*brn$time
brn$Status <- rep(0, dim(brn)[1])
brn$Status[!oo] <- 1

brn <- brn %>% select(-c(Patient.ID, Relapse.Free.Status, Relapse.Free.Status..Months., Overall.Survival.Status, Overall.Survival..Months., TMB..nonsynonymous.)) 

####OUR DATA IS brn






#############baseline characteristics of the data
brn %>% group_by(Status) %>% summarise(avg = mean(age))
brn %>% group_by(Status) %>% summarise(avg = sd(age))
brn %>% group_by(Status) %>% summarise(avg = mean(Lymph.nodes.examined.positive))
brn %>% group_by(Status) %>% summarise(avg = sd(Lymph.nodes.examined.positive))
brn %>% group_by(Status) %>% summarise(avg1= mean(Mutation.Count), avg = sd(Mutation.Count))
brn %>% group_by(Status) %>% summarise(avg1= mean(Nottingham.prognostic.index), avg = sd(Nottingham.prognostic.index))
brn %>% group_by(Status) %>% summarise(avg1= mean(Tumor.Size), avg = sd(Tumor.Size))
brn %>% group_by(Status) %>% summarise(percent = Radio.Therapy/sum(Radio.Therapy))
table(brn$ER.Status, brn$Status)

##------------------------------------------------------------
## compare RF-SRC to Cox regression
## Illustrates C-index and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------

library(survival)
library(pec)
library(party)
library(randomForestSRC)
library(prodlim)
library(ipred)
library(vip)



  trainrows <-  sample(1:nrow(brn), .7*nrow(brn))
  
  # GBSG2 <- GBSG2 %>% mutate(event = factor(cens, 
  #                                         labels =  c("censored", "failed")))
  #table(GBSG2$event)
  #GBSG2 <- GBSG2[,-10]
  surv.f.br <- as.formula(Surv(time, Status) ~.)
  pec.f.br <- as.formula(Hist(time, Status) ~ 1)
  
  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj.br <- coxph(surv.f.br, data = brn, x = TRUE)
  rfsrc.obj.br <- rfsrc(surv.f.br, brn, ntree = 1000)
  fitcforest.br <- pecCforest(surv.f.br, data = brn, controls = cforest_classical(ntree = 1000))
  ### estimate conditional Kaplan-Meier curves
  bst.br <- cforest(Surv(time, Status) ~ ., data = brn,
                 control = cforest_unbiased(ntree = 50))
  
  treeresponse(bst, newdata = GBSG2[1:2,], OOB = TRUE)
  ### if you can't resist to look at individual trees ...
  party:::prettytree(bst@ensemble[[1]], names(bst@data@get("input")))
  vi <- varimp(bst.br)
  ### compare variable importances and absolute z-statistics
  layout(matrix(1:2))
  barplot(vi)
  barplot(abs(summary(coxph(Surv(time, cens) ~ ., data = GBSG2))$coeff[,"z"])) ### looks more or less the same
  
  
  
  plot.variable(rfsrc.obj, surv.type = "rel.freq", partial = TRUE)
  vip(bst.br)#Variable importance for survival forests (cforest)
  
  bst.importance <- varimp(bst.br, conditional = TRUE)
  names(bst.importance) <- c("Cell", "Chemo", "ER", "Neo", "HER2", "HormoneT", "Menostate","Lymph",
                         "Mut.count","Notti.Prog", "PR", "Radio", "Tsize", "Tgrade", "Age")
  
  par(mfrow=c(1,1))
  dotchart(sort(bst.importance), main = "Conditional Importance of Variables for METABRIC Data")
  
  
  rsf.cindex.br <- get.cindex(time = brn$time, censoring = brn$Status, predicted = rfsrc.obj.br$predicted.oob)
  
  calPlot(rfsrc.obj, data = )
  
  newData.br <- data.frame(Cellularity=factor("Moderate",levels=c("Low", "Moderate", "High")),
                        Chemotherapy=factor("NO",levels=c("NO", "YES")),
                        ER.Status=factor("Negative",levels=c("Negative", "Positive")),
                        Neoplasm.Histologic.Grade=factor(2, levels=c(1,2,3)),
                        HER2.Status=factor("Negative",levels=c("Negative", "Positive")),
                        Hormone.Therapy=factor("NO",levels=c("NO", "YES")),
                        age=as.integer(c(24,53,73)),
                        Inferred.Menopausal.State=factor(cbind("Pre", "Post", "Post"),levels=c("Post","Pre")),
                        Lymph.nodes.examined.positive =rep(3, 3),
                        Mutation.Count=rep(4,3),
                        Nottingham.prognostic.index=rep(4.12,3),
                        PR.Status=factor("Negative",levels=c("Negative", "Positive")),
                        Radio.Therapy=factor("NO",levels=c("NO", "YES")),
                        TMB..nonsynonymous. = rep(7.16, 3),
                        Tumor.Size=as.numeric(rep(25, 3)),
                        Tumor.Stage=factor(2,levels=c(1,2,3)))
  print(newData.br)
  
  
  
  
  cox.obj2.br <- coxph(surv.f.br, data = brn[trainrows,], x = TRUE)
  rfsrc.obj2.br <- rfsrc(surv.f.br, brn[trainrows,], ntree = 1000)
  fitcforest2.br <- pecCforest(surv.f.br, data = brn[trainrows,], controls = cforest_classical(ntree = 1000))
  
  
  
  
  
  
  cforestImpPlot <- function(x) {
    cforest_importance <<- v <- varimp(x)
    dotchart(v[order(v)])
  }
  cforestImpPlot(bst.br)
  
  ########################################################################################
  #################predicted survival curves for newdata#############################
  ########################################################################################
  
  pcox <- predictSurvProb(cox.obj.br, newdata = newData.br, times = quantile(brn$time))
  prsf <- predictSurvProb(rfsrc.obj.br, newdata = newData.br, times = quantile(brn$time))
  pcf <- predictSurvProb(fitcforest.br, newdata = newData.br, times = quantile(brn$time))
  psp <- round(cbind(pcox,prsf, pcf),2)
  newNames <- paste("newData", 1:3)
  mat <- cbind( newData$age, psp)
  colnames(mat) <- c( "Age", rep("Coxph", 5), rep("rsf",5), rep("cforest",5))
  mat1 <- data.frame(newNames, mat)
  
  print(mat1)
  
  
  par(mfrow = c(1, 3))
  lapply(1:3, function(x) {
    plotPredictSurvProb(cox.obj.br, newdata = newData.br[x, ], col = 1)
    plotPredictSurvProb(rfsrc.obj.br, newdata = newData.br[x, ], add = TRUE, col = 2)
    plotPredictSurvProb(fitcforest.br, newdata = newData.br[x, ], add = TRUE,col = 3)
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
  set.seed(13)
  ptm <- proc.time()
  prederror.brn <- pec(list(cox.obj2.br,rfsrc.obj2.br, fitcforest2.br), data = brn[-trainrows,], formula = pec.f.br,
                         splitMethod = "bootcv", B = 500, M=200)
  print(prederror.brn, times = seq(8, 2659,50))
  plot(prederror.brn)
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  
  
  set.seed(13)
  ptm <- proc.time()
  prederror.brn.632 <- pec(list(cox.obj2.br,rfsrc.obj2.br, fitcforest2.br), data = brn[-trainrows,], formula = pec.f.br,
                             splitMethod = "Boot632", B = 500, M=200)
  print(prederror.brn.632, times = seq(8, 2659,50))
  plot(prederror.brn.632)
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  
  
  plot(prederror.brn, legend=FALSE)
  legend("bottomright", legend=c("Reference", "Cox PH","RSF", "Cforest"),
         col=c("black", "red","green", "blue"), lty=1)
  title(main = "METABRIC bootstrap estimates")
  plot(prederror.brn.632, legend=FALSE)
  legend("bottomright", legend=c("Reference", "Cox PH","RSF", "Cforest"),
         col=c("black", "red","green", "blue"), lty=1)
  title(main = "METABRIC bootstrap.632 estimates")
  ########################################################################################
  #########################################################################################
  ########################################################################################
  
  
  
  
  #########################################################################################
  ###########cindex##############################################################################
  ########################################################################################
  
  # compute the bootstrap-crossvalidation estimate of
  # the C-index at different time points #
  ptm <- proc.time()
  set.seed(13)
  bcvCindex.br <- pec::cindex(list("Cox"=cox.obj.br,
                                "RSF"=rfsrc.obj.br,
                                "Cforest"= fitcforest.br), formula=surv.f.br,
                           data=brn,
                           splitMethod="Bootcv",
                           B = 500, M=200, eval.times=c(50,quantile(brn$time)[-1]))
  bcvCindex.632.br <- pec::cindex(list("Cox"=cox.obj.br,
                                    "RSF"=rfsrc.obj.br,
                                    "Cforest"= fitcforest.br), formula=surv.f.br,
                               data=brn,
                               splitMethod="Boot632",
                               B = 500, M=200, eval.times=c(50,quantile(brn$time)[-1]))
  print(bcvCindex.br)
  par(mfrow = c(2,2))
  
  plot(bcvCindex.br, legend=FALSE)
  legend("topleft", legend=c("Cox PH","RSF", "Cforest"),
         col=c("black", "red","green"), lty=1)
  title(main = "METABRIC bootstrap estimates")
  
  print(bcvCindex.632.br)
  plot(bcvCindex.632.br, legend=FALSE)
  legend("topleft", legend=c("Cox PH","RSF", "Cforest"),
         col=c("black", "red","green"), lty=1)
  title(main = "METABRIC bootstrap .632 estimates")
  stoptime<- proc.time() - ptm # Stop the clock
  print(stoptime[3]/60)
  ########################################################################################
  #########################################################################################
  ########################################################################################
  
  
  
  
  
  #Rkare <- R2(prederror.GBSG2,times=seq(8, 2659,50))
  
  ######calibration plot cizdirmeye calistim ama olmadi....
  #veriyi train diye bolup oyle model kurdum hata bulmaya calisiyorum
  set.seed(14)
  prederror2lam <- pec(list(cox.obj2,rfsrc.obj2, fitcforest2), data = GBSG2[trainrows,], formula = pec.f,
                    splitMethod = "bootcv", B = 50)
  print(prederror2)
  plot(prederror2)
  

  
  
  # oo <- subsample(rfsrc.obj)
  # # take a delete-d-jackknife procedure for example
  # vimpCI <- extract.subsample(oo)$var.jk.sel.Z
  # # Confidence Intervals for VIMP
  # plot.subsample(oo)
  # # take the variable "Month" for example for partial plot
  # plot.variable(rfsrc.obj, xvar.names = "days", partial = TRUE)
  
  
  
  
  
  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj.brr <- rfsrc(surv.f.br, brn)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(brn), nrow(brn), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f.br, brn[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj.br)) {
      get.cindex(brn$time[-train], brn$Status[-train], predict(cox.obj.br, brn[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj.br$err.rate[rfsrc.obj.br$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")







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


brn.new <- brn
colnames(brn.new) <- c("Cell", "Chemo", "ER", "Neo", "HER2", "HormoneT", "Menostate","Lymph",
                   "Mut.count","Notti.Prog", "PR", "Radio", "Tsize", "Tgrade", "Age", "time", "Status")

brn.log1 <- rfsrc(Surv(time, Status) ~ ., brn.new, 
                   importance="permute")

brn.cr <- rfsrc(Surv(time, Status) ~ ., brn.new)


plot(brn.log1)# bunu ekledim sonuclara



## log-rank cause-2 specific splitting and targeted VIMP for cause 2
#GBSG2.log2 <- rfsrc(Surv(time, cens) ~ ., GBSG2, 
#                 splitrule = "logrank", cause = c(0,1), importance = TRUE)

## extract VIMP from the log-rank forests: event-specific
## extract minimal depth from the Gray log-rank forest: non-event specific
var.perf.br <- data.frame(md = max.subtree(brn.cr)$order[, 1],
                       vimp = 100 * brn.log1$importance)
lap.br <- print(var.perf.br[order(var.perf.br$vimp), ])
stargazer(cox.obj,summary=TRUE, type = "text", title = "var per", out="cox.obj.txt")


plot(GBSG2.log1)


p<-ggplot(data=lap.br, aes(x=md, y=vimp)) + geom_bar(stat="identity")+  geom_text(aes(label=vimp), vjust=1.6, color="white", size=3.5)+
  theme_minimal()
p + coord_flip()


# Inside bars
ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()



############print results
A <- summary(cox.obj.br)
write.table(round(A$coefficients,4), file = "cox.obj.br.txt", sep = ",", quote = FALSE)

B <- as.data.frame(t(mat1))
write.table(B, "metabric.newdata.txt", sep = ",", quote = FALSE)

A2 <- summary(cox.obj2.br)
write.table(round(A2$coefficients,4), file = "cox.obj.br2.txt", sep = ",", quote = FALSE)


save.image("cforest.rf_br.RData")

save.image("cforest.rf_brandGBSG2.RData")
