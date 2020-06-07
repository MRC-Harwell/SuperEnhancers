#!/bin/Rscript
library(randomForest)
library(ROCR)
library(caret)
library(e1071)
library(PRROC)
library(doMC)
registerDoMC(cores = 10)

# df1 = ML table matrix with data (TSRE;PPI;TF;EXP)
mp = "Mgi_phenotypes/MATRIX.txt" # phenotype matrix to be used as response variable (binary : 0/1)

df1 = cbind(df0,df3,df4,df5,df6,df7,df8)

i = 2 # column number of phenotype to predict in phenotype matrix
term = colnames(df2[i])
print(term)
df1$MP = df2[,i]
df1$MP = ifelse(df1$MP == 0, "absent", df1$MP)
df1$MP = ifelse(df1$MP == 1, "present", df1$MP)


cbind.fill <- function(...){
	nm <- list(...)
	nm <- lapply(nm, as.matrix)
	n <- max(sapply(nm, nrow))
	do.call(cbind, lapply(nm, function (x)
			rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


cm = data.frame()
auc.roc = data.frame()
auc.pr = data.frame()
dec.accu = data.frame()
dec.gini = data.frame()
res.pred = list()
res.labels = list()

counter = 0
for(r in 1:10){
	# set seed here as well
	flds <- createFolds(factor(df1$MP), k = 5, list = TRUE, returnTrain = FALSE)
	for(j in 1:5){
		counter = counter + 1
		test = df1[unlist(flds[j]),]
		train = df1[-unlist(flds[j]),]

		track = paste("REPEAT: ",r," -> Fold number: ",j,sep="")
		cat("\n*************************************\n")
		print(track)
		cat("***************************************\n")

		set.seed(10) #run 1
		seeds <- vector(mode = "list", length = 26)
		m = round(sqrt(ncol(train)))
		for(i in 1:25) seeds[[i]]<- sample.int(n=1000, m)
		seeds[[26]]<-sample.int(1000, 1)

		cntrl = trainControl(method="repeatedcv", number=5, repeats=5,summaryFunction = twoClassSummary,sampling = "down", allowParallel = TRUE, classProbs = TRUE, seeds = seeds)
		rf_model<-train(factor(MP)~.,data=train,method="rf",trControl=cntrl, metric = "ROC",importance = TRUE, preProc = c("center","scale"))

		############### for predicting training ################
		# train$predicted.response = predict(rf_model, train)
		# confuseMat.train = confusionMatrix(data=train$predicted.response, reference=train$MP, positive = "present")
		# confuseMat.train
		# varImp(rf_model, scale=F)

		################ importance of variables #############
		var.imp = importance(rf_model$finalModel)
		dec.accu = cbind.fill(dec.accu,var.imp[,3])
		dec.gini = cbind.fill(dec.gini,var.imp[,4])

		################ predicting test data and producing confusion Matrix #######################
		test$predicted.response <- predict(rf_model, test)
		confuseMat.test = confusionMatrix(data=test$predicted.response,reference=test$MP, positive = "present")
		results = t(data.frame(cbind(t(confuseMat.test$byClass),t(confuseMat.test$overall))))
		accuracy = round(results[12],2)
		sens = round(results[1],2)
		spec = round(results[2],2)
		table = confuseMat.test$table

		cm = cbind.fill(cm,results[,1])

		###################### validating test data ####################

		## For producing probabilities ####
		test.pr = predict(rf_model, test, type = "prob")[,2]
		test.pred = prediction(test.pr,test$MP)

		test.AUC=performance(test.pred,"auc") #Calculate the AUC value
		AUC=test.AUC@y.values[[1]]

		# storing ROCR curve values
		auc.roc[counter,"ROCAUC"] = AUC
		res.pred[[counter]] <- test.pr
		res.labels[[counter]] <- test$MP

		#PR curve
		fg = test.pr[test$MP == "present"]
		bg = test.pr[test$MP == "absent"]
		pr = pr.curve(scores.class0 = fg, scores.class1 = bg, curve=T)
		prauc = pr$auc.integral

		# storing PR curve values
		auc.pr[counter,"PRAUC"] = prauc
	}
}
# writing the results to files
# Mean decrease accuracy
accu.df <- data.frame(unlist(dec.accu),stringsAsFactors=FALSE)
accu.df.mean = rowMeans(accu.df)
accu.df = transform(accu.df, sd=apply(accu.df,1, sd))
accu.df$mean = accu.df.mean
file8 = paste(term,".decAccuracy.txt",sep="")
write.table(accu.df,file = file8, sep="\t",quote=F)

# Mean decrease Gini
gini.df <- data.frame(unlist(dec.gini),stringsAsFactors=FALSE)
gini.df.mean = rowMeans(gini.df)
gini.df = transform(gini.df, sd=apply(gini.df,1, sd))
gini.df$mean = gini.df.mean
file9 = paste(term,".decGini.txt",sep="")
write.table(gini.df,file = file9, sep="\t",quote=F)

# stats from confusion matrix
cm.df <- data.frame(unlist(cm),stringsAsFactors=FALSE)
cm.df.mean = rowMeans(cm.df)
cm.df = transform(cm.df, sd=apply(cm.df,1, sd))
cm.df$mean = cm.df.mean
file3 = paste(term,".results.txt",sep="")
write.table(cm.df,file = file3, sep="\t",quote=F,row.names=F)

# auc for ROC and PR
file4 = paste(term,".auc.txt",sep="")
write.table(auc.roc,file=file4,sep="\t",quote=F,row.names=F, col.names=T)
file5 = paste(term,".prauc.txt",sep="")
write.table(auc.pr,file=file5,sep="\t",quote=F,row.names=F, col.names=T)

########################################################################
######################## ROCR curve ####################################
pred.all = prediction(res.pred,res.labels)
perf.all = performance(pred.all,"tpr","fpr")

# function to average results, made by modifying function from ROCR package
avg.results <- function (perf)
{
	## for infinite cutoff, assign maximal finite cutoff + mean difference
	  ## between adjacent cutoff pairs
	  if (length(perf@alpha.values)!=0) perf@alpha.values <-
	    lapply(perf@alpha.values,
	           function(x) { isfin <- is.finite(x);
	                         x[is.infinite(x)] <-
	                           (max(x[isfin]) +
	                            mean(abs(x[isfin][-1] -
	                                     x[isfin][-length(x[isfin])])));
	                         x } )
	  ## remove samples with x or y not finite
	  for (i in 1:length(perf@x.values)) {
	      ind.bool <- (is.finite(perf@x.values[[i]]) &
	                   is.finite(perf@y.values[[i]]))

	      if (length(perf@alpha.values)>0)
	        perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]

	      perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
	      perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
	  }

    perf.sampled <- perf
    alpha.values <- rev(seq(min(unlist(perf@alpha.values)), max(unlist(perf@alpha.values)),
        length = max(sapply(perf@alpha.values, length))))
    for (i in 1:length(perf.sampled@y.values)) {
        perf.sampled@x.values[[i]] <- approxfun(perf@alpha.values[[i]],
            perf@x.values[[i]], rule = 2, ties = mean)(alpha.values)
        perf.sampled@y.values[[i]] <- approxfun(perf@alpha.values[[i]],
            perf@y.values[[i]], rule = 2, ties = mean)(alpha.values)
    }
    perf.avg <- perf.sampled
	perf.avg.data <<- perf.avg
    perf.avg@x.values <- list(rowMeans(data.frame(perf.avg@x.values)))
    perf.avg@y.values <- list(rowMeans(data.frame(perf.avg@y.values)))
    perf.avg@alpha.values <- list(alpha.values)
	perf.rocr.avg <<- perf.avg
}

avg.results(perf.all)

rocs.x = data.frame()
for(q in 1:length(perf.avg.data@x.values)){
	values = as.data.frame(unlist(perf.avg.data@x.values[q]))
	colnames(values) = c("x")
	rocs.x = rbind(rocs.x,values)
}

rocs.y = data.frame()
for(q in 1:length(perf.avg.data@y.values)){
	values = as.data.frame(unlist(perf.avg.data@y.values[q]))
	type1 = paste("fold",q, sep="")
	type2 = "cv"
	values$type1 = type1
	values$type2 = type2
	colnames(values) = c("y","type1","type2")
	rocs.y = rbind(rocs.y,values)
}

roc.all = cbind(rocs.x,rocs.y)
file6 = paste(term,".rocValues.cv.txt",sep="")
write.table(roc.all,file=file6,sep="\t",quote=F,row.names=F)

rocs.y.wide = data.frame(perf.avg.data@y.values)
avg.y = rowMeans(rocs.y.wide)
rocs.y.wide = transform(rocs.y.wide, sd=apply(rocs.y.wide,1, sd))
rocs.y.wide$mean = avg.y
file12 = paste(term,".rocMeanSd.txt",sep="")
write.table(rocs.y.wide,file=file12,sep="\t",quote=F,row.names=F)

mean.x = as.data.frame(unlist(perf.rocr.avg@x.values))
colnames(mean.x) = c("x")
mean.y = as.data.frame(unlist(perf.rocr.avg@y.values))
mean.y$type1 = "mean"
mean.y$type2 = "mean"
colnames(mean.y) = c("y","type1","type2")
mean.all = cbind(mean.x,mean.y)
file7 = paste(term,".rocValues.mean.txt",sep="")
write.table(mean.all,file=file7,sep="\t",quote=F,row.names=F)

roc.all = rbind(roc.all,mean.all)
m.rauc = round(mean(auc.roc$ROCAUC),3)
sd.rauc = round(sd(auc.roc$ROCAUC),3)

a1 = paste("mean"," (AUC=",m.rauc,"\u00B1",sd.rauc,")", sep="")


######################################################################
######################## PR curve ####################################
perf.pr.all = performance(pred.all,"prec","rec")
avg.results(perf.pr.all)

pr.x = data.frame()
for(q in 1:length(perf.avg.data@x.values)){
	values = as.data.frame(unlist(perf.avg.data@x.values[q]))
	colnames(values) = c("x")
	pr.x = rbind(pr.x,values)
}

pr.y = data.frame()
for(q in 1:length(perf.avg.data@y.values)){
	values = as.data.frame(unlist(perf.avg.data@y.values[q]))
	type1 = paste("fold",q, sep="")
	type2 = "cv"
	values$type1 = type1
	values$type2 = type2
	colnames(values) = c("y","type1","type2")
	pr.y = rbind(pr.y,values)
}

pr.all = cbind(pr.x,pr.y)
file10 = paste(term,".prValues.cv.txt",sep="")
write.table(pr.all,file=file10,sep="\t",quote=F,row.names=F)

pr.y.wide = data.frame(perf.avg.data@y.values)
avg.pr.y = rowMeans(pr.y.wide)
pr.y.wide = transform(pr.y.wide, sd=apply(pr.y.wide,1, sd))
pr.y.wide$mean = avg.pr.y
file12 = paste(term,".prMeanSd.txt",sep="")
write.table(pr.y.wide,file=file12,sep="\t",quote=F,row.names=F)

mean.pr.x = as.data.frame(unlist(perf.rocr.avg@x.values))
colnames(mean.pr.x) = c("x")
mean.pr.y = as.data.frame(unlist(perf.rocr.avg@y.values))
mean.pr.y$type1 = "mean"
mean.pr.y$type2 = "mean"
colnames(mean.pr.y) = c("y","type1","type2")
mean.pr.all = cbind(mean.pr.x,mean.pr.y)
file11 = paste(term,".prValues.mean.txt",sep="")
write.table(mean.pr.all,file=file11,sep="\t",quote=F,row.names=F)

pr.all = rbind(pr.all,mean.pr.all)
m.prauc = round(mean(auc.pr$PRAUC),3)
sd.prauc = round(sd(auc.pr$PRAUC),3)

a2 = paste("mean"," (AUC=",m.prauc,"\u00B1",sd.prauc,")", sep="")

#########################################
