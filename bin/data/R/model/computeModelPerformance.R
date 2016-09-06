library("ROCR")
args <- commandArgs(TRUE)
inFile<-args[1]
outFile<-args[2]
expThreshold<-args[3]

R_Model <- read.table(inFile, as.is= TRUE, sep = "\t", header=TRUE)

#FABIO TODO: use expThreshold to compute categories 

table<-R_Model[complete.cases(R_Model),]
fit = lm(table$experimental~table$prediction)
m = summary(fit)
SE <- m$sigma
R2 <- m$r.squared
Tot <- (nrow(table)-1)
Pos <- length(which(table$experimental>=expThreshold))
Neg <- length(which(table$experimental<expThreshold))
Class<-findInterval(table$experimental,expThreshold)
pred<-prediction(table$prediction, Class)
perf<-performance(pred,"auc")
AUC = perf@y.values
perf<-performance(pred,"tpr","tnr")
listTPR <-simplify2array(perf@y.values)
listTNR <-simplify2array(perf@x.values)
perf<-performance(pred,"acc")
listACC <- simplify2array(perf@y.values)
AvgAccuracy=data.frame(listTPR, listTNR)
AvgAccuracy2=data.frame(AvgAccuracy, listACC)
AvgAccuracyValue = max(rowMeans(AvgAccuracy, na.rm = TRUE))
bestindex = which.max(rowMeans(AvgAccuracy, na.rm = TRUE))
listCutOff <-simplify2array(perf@x.values)
TPR<-listTPR[bestindex]
TNR<-listTNR[bestindex]
CCM<-rowMeans(AvgAccuracy, na.rm = TRUE)[bestindex]
CutOff<-listCutOff[bestindex]
perf<-performance(pred,"ppv","npv")
listPPV <-simplify2array(perf@y.values)
listNPV <-simplify2array(perf@x.values)
PPV<-listPPV[bestindex]
NPV<-listNPV[bestindex]
Value <- c(signif(simplify2array(AUC),3),signif(simplify2array(TPR),3),signif(simplify2array(TNR),3),signif(simplify2array(CCM),3),signif(simplify2array(PPV),3),signif(simplify2array(NPV),3),signif(simplify2array(CutOff),3),signif(simplify2array(SE),3),signif(simplify2array(R2),3),simplify2array(Tot),simplify2array(Pos),simplify2array(Neg))
Title <- c("AUC","TPR","TNR","CCM","PPV","NPV","CutOff","SE","R2","Tot","Pos","Neg")
Output = data.frame(Title, Value)

write.table(Output, outFile ,sep = "\t") 

quit("no")
