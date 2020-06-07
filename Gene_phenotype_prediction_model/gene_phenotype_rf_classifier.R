#!/bin/Rscript
library(randomForest)
library(ROCR)
library(caret)
library(e1071)
library(PRROC)
library(doMC)
registerDoMC(cores = 13)
library(ggplot2)

#qsub -cwd -j y -b yes -N mp.5367 -P NGS -o /NGS/users/Sid/ENHANCERS/SuperEnh_27ac/Great_targets/Strict_targets/Final_predictionModels_3/Rf_exp_tsre_tsrePpi_ppi_tf/Logs/ -q big.q@kansas.mrch.har.mrc.ac.uk -pe big 16 /R/R-IMPC-sid/bin/R CMD BATCH /NGS/users/Sid/SCRIPTS/Chromatin_segmentation/RandomForest/randomFr_working_singleMp_v2.R MP.0005367.Rout

enhancer = "/NGS/users/Sid/ENHANCERS/SuperEnh/AllGenes_cumStrongEnh/MATRIX.txt"
promoter =  "/NGS/users/Sid/PROMOTERS/AllGenes_cumStrongProm/MATRIX.txt"
mp = "/NGS/users/Sid/ENHANCERS/SuperEnh/Mgi_phenotypes/MATRIX.txt"
expression = "/NGS/users/Sid/ENHANCERS/SuperEnh/AllGenes_exp/MATRIX.txt"
enh_ppi = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/AllGenes_ppiEnh/MATRIX.txt"
prom_ppi =  "/NGS/users/Sid/PROMOTERS/AllGenes_ppiProm/MATRIX.txt"
mp_ppi = "/NGS/users/Sid/ENHANCERS/SuperEnh/AllGenes_ppiMP/MP:0005397.txt"
tf = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/TF_features"

df0 = read.table(expression,header=T,sep="\t")
rownames(df0) = df0$Gene
colnames(df0) = c("Gene","BAT_Exp","BmarrowDm_Exp","Bmarrow_Exp","Cerebellum_Exp","CH12_Exp","Cortex_Exp","Esb4_Exp","Es.E14_Exp","Heart_Exp","Kidney_Exp","Limb_Exp","Liver_Exp","Lung_Exp","MEF_Exp","MEL_Exp","OlfactoryBulb_Exp","Placenta_Exp","SmallIntestine_Exp","Spleen_Exp","Testis_Exp","Thymus_Exp","Wbrain_Exp")
df0 = df0[,-1]


##########data = subset(df1,select = c("BAT","Bmarrow","Heart","Kidney","Testis","Wbrain"))

df2 = read.table(mp,header=T,sep="\t")

df3 = read.table(enhancer,header=T,sep="\t")
colnames(df3) = c("Gene","BAT_Enh","BmarrowDm_Enh","Bmarrow_Enh","Cerebellum_Enh","CH12_Enh","Cortex_Enh","Esb4_Enh","Es.E14_Enh","Heart_Enh","Kidney_Enh","Limb_Enh","Liver_Enh","Lung_Enh","MEF_Enh","MEL_Enh","OlfactoryBulb_Enh","Placenta_Enh","SmallIntestine_Enh","Spleen_Enh","Testis_Enh","Thymus_Enh","Wbrain_Enh")
rownames(df3) = df3$Gene
df3 = df3[,-1]
#df3[] <- lapply(df3, as.factor)

df4 = read.table(promoter,header=T,sep="\t")
colnames(df4) = c("Gene","BAT_Prom","BmarrowDm_Prom","Bmarrow_Prom","Cerebellum_Prom","CH12_Prom","Cortex_Prom","Esb4_Prom","Es.E14_Prom","Heart_Prom","Kidney_Prom","Limb_Prom","Liver_Prom","Lung_Prom","MEF_Prom","MEL_Prom","OlfactoryBulb_Prom","Placenta_Prom","SmallIntestine_Prom","Spleen_Prom","Testis_Prom","Thymus_Prom","Wbrain_Prom")
rownames(df4) = df4$Gene
df4 = df4[,-1]
#df4[] <- lapply(df4, as.factor)

df5 = read.table(enh_ppi,header=T,sep="\t")
colnames(df5) = c("Gene","BAT_EnhPPI","BmarrowDm_EnhPPI","Bmarrow_EnhPPI","Cerebellum_EnhPPI","CH12_EnhPPI","Cortex_EnhPPI","Esb4_EnhPPI","Es.E14_EnhPPI","Heart_EnhPPI","Kidney_EnhPPI","Limb_EnhPPI","Liver_EnhPPI","Lung_EnhPPI","MEF_EnhPPI","MEL_EnhPPI","OlfactoryBulb_EnhPPI","Placenta_EnhPPI","SmallIntestine_EnhPPI","Spleen_EnhPPI","Testis_EnhPPI","Thymus_EnhPPI","Wbrain_EnhPPI")
rownames(df5) = df5$Gene
df5 = df5[,-1]

df6 = read.table(prom_ppi,header=T,sep="\t")
colnames(df6) = c("Gene","BAT_PromPPI","BmarrowDm_PromPPI","Bmarrow_PromPPI","Cerebellum_PromPPI","CH12_PromPPI","Cortex_PromPPI","Esb4_PromPPI","Es.E14_PromPPI","Heart_PromPPI","Kidney_PromPPI","Limb_PromPPI","Liver_PromPPI","Lung_PromPPI","MEF_PromPPI","MEL_PromPPI","OlfactoryBulb_PromPPI","Placenta_PromPPI","SmallIntestine_PromPPI","Spleen_PromPPI","Testis_PromPPI","Thymus_PromPPI","Wbrain_PromPPI")
rownames(df6) = df6$Gene
df6 = df6[,-1]

df7 = read.table(mp_ppi,header=F,sep="\t")
colnames(df7) = c("Gene","PPI_score")
rownames(df7) = df7$Gene
df7 = df7[,-1, drop=F]

# TF ##
filenames <- list.files(tf, pattern="*.txt", full.names=TRUE)
tf_path = "/NGS/users/Sid/ENHANCERS/SuperEnh_27ac/TF_features/"
readIn <- function(x) {
  data.f = read.table(x,header=F,sep="\t", stringsAsFactor=F)
  data.f = data.f[,2,drop=F]
  name = sub(tf_path,"",x)
  name = sub(".txt","",name)
  colnames(data.f) = c(name)
  return(data.f)
}

head = head(filenames)
ldf.tmp <- lapply(head, function(x) read.table(x,header=F,sep="\t", stringsAsFactor=F))
tmp.df = ldf.tmp[[1]]

ldf.tf <- lapply(filenames, readIn)
df8 <- do.call(cbind, ldf.tf)
rownames(df8) = tmp.df$V1


df1 = cbind(df0,df3,df4,df5,df6,df7,df8)  # all + tf
#df1 = cbind(df0,df3,df4,df5,df6,df7)  # all
#df1 = cbind(df0,df7) # exp+ppi
#df1 = df8  #tf
#df1 = cbind(df3,df4,df8)  # tsre + tf
#df1 = df0  #exp

#################################################################################################

i = 2
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
		#confuseMat.test
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

# function to average results, made by stealing and hijacking function from ROCR package
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

# Plotting ROCR curve
# file1 = paste(term,".roc.png",sep="")
# #title = paste("Accuracy=",accuracy,", Sensitivity=",sens,", Specificity=",spec,sep="")
# png(file1,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
# ggplot() +
# geom_line(data = roc.all[roc.all$type2 == "cv",], aes(x=x, y=y, group = type1, color = type2), size = 0.2, alpha = 0.3) +
# geom_line(data = roc.all[roc.all$type2 == "mean",], aes(x=x, y=y, color = type2), size = 0.5) +
# geom_abline(intercept = 0, linetype = "dashed", color = "black", alpha = 0.3)+
# theme_bw()+
# scale_colour_manual(
# 	values= c("grey","blue"),
# 	breaks = c("cv","mean"),
# 	labels = c("cross validation",a1))+
# geom_ribbon(data = rocs.y.wide, aes(roc.all[roc.all$type2 == "mean",]$x,ymax = mean + sd, ymin = mean - sd), alpha = 0.2, fill = "blue")+
# theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
# theme(plot.title = element_text(size=11, hjust = 0.5),
# 	legend.key.size = unit(0.30,"cm"),
# 	legend.title = element_blank(),
# 	legend.text = element_text(size = 7),
# 	legend.position = c(0.80,0.08),
# 	legend.background = element_rect(fill = "transparent",colour = "NA"),
# 	axis.text.x = element_text(size=8),
# 	axis.text.y = element_text(size=8),
# 	panel.border = element_rect(colour="BLACK",size=0.4),
# 	axis.title.x = element_text(size=11, vjust = 0.1),
# 	axis.title.y = element_text(size=11,angle = 90, vjust = 0.1),
# 	panel.background = element_rect(fill="transparent"),
# 	plot.background = element_rect(fill = "transparent",colour = NA)
# ) +
# scale_y_continuous(name="True positive rate") +
# scale_x_continuous(name = "False positive rate")
# #ggtitle(title)
# dev.off()

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

# Plotting PR curve
# file2 = paste(term,".pr.png",sep="")
# #title = paste("Accuracy=",accuracy,", Sensitivity=",sens,", Specificity=",spec,sep="")
# png(file2,bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
# ggplot() +
# geom_line(data = pr.all[pr.all$type2 == "cv",], aes(x=x, y=y, group = type1, color = type2), size = 0.2, alpha = 0.3) +
# geom_line(data = pr.all[pr.all$type2 == "mean",], aes(x=x, y=y, color = type2), size = 0.5) +
# theme_bw()+
# scale_colour_manual(
# 	values= c("grey","blue"),
# 	breaks = c("cv","mean"),
# 	labels = c("cross validation",a2))+
# geom_ribbon(data = pr.y.wide, aes(pr.all[pr.all$type2 == "mean",]$x,ymax = mean + sd, ymin = mean - sd), alpha = 0.2, fill = "blue")+
# theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
# theme(plot.title = element_text(size=11, hjust = 0.5),
# 	legend.key.size = unit(0.30,"cm"),
# 	legend.title = element_blank(),
# 	legend.text = element_text(size = 7),
# 	legend.position = c(0.80,0.96),
# 	legend.background = element_rect(fill = "transparent",colour = "NA"),
# 	axis.text.x = element_text(size=8),
# 	axis.text.y = element_text(size=8),
# 	panel.border = element_rect(colour="BLACK",size=0.4),
# 	axis.title.x = element_text(size=11, vjust = 0.1),
# 	axis.title.y = element_text(size=11,angle = 90, vjust = 0.1),
# 	panel.background = element_rect(fill="transparent"),
# 	plot.background = element_rect(fill = "transparent",colour = NA)
# ) +
# scale_y_continuous(name="Precision") +
# scale_x_continuous(name = "Recall")
# #ggtitle(title)
# dev.off()

################ Hacking the ROCR package ###########
###### trace(.performance.plot.threshold.avg, edit=TRUE)
######  line added: perf.rocr.avg <<- perf.avg ###




	# ##############################################################
	# ##############################################################
	# ######### For exploring False positives ######################
	# ##############################################################
	# set.seed(10) #run 1
	# seeds <- vector(mode = "list", length = 26)
	# m = round(sqrt(ncol(df1)))
	# for(i in 1:25) seeds[[i]]<- sample.int(n=1000, m)
	# seeds[[26]]<-sample.int(1000, 1)
	#
	# cntrl = trainControl(method="repeatedcv", number=5, repeats=5,summaryFunction = twoClassSummary,sampling = "down", allowParallel = TRUE, classProbs = TRUE, seeds = seeds)
	# rf_model<-train(factor(MP)~.,data=df1,method="rf",trControl=cntrl, metric = "ROC",importance = TRUE, preProc = c("center","scale"))
	#
	# ################ predicting using ALL data and producing confusion Matrix #######################
	# df1$response <- predict(rf_model, df1)
	# confuseMat.test = confusionMatrix(data=df1$response,reference=df1$MP, positive = "present")
	# #confuseMat.test
	# file11 = paste(term,".ALL.results.txt",sep="")
	# file12 = paste(term,".ALL.table.txt",sep="")
	# results = t(data.frame(cbind(t(confuseMat.test$byClass),t(confuseMat.test$overall))))
	# accuracy = round(results[12],2)
	# sens = round(results[1],2)
	# spec = round(results[2],2)
	# table = confuseMat.test$table
	# write.table(results,file=file11,sep="\t",quote=F,col.names=F)
	# write.table(table,file=file12,sep="\t",quote=F)
	#
	# df1$prob = predict(rf_model, df1, type = "prob")[,2]
	# predictions = subset(df1,select = c("MP","response","prob"))
	# predictions$gene = rownames(predictions)
	# all.p = predictions[predictions$response == "present",]
	# fp = all.p[all.p$MP == "absent",]
	# fp.sort = fp[order(fp$prob, decreasing=T),]
	#
	# file13 = paste(term,".ALL.FP.txt",sep="")
	# write.table(fp.sort,file=file13,sep="\t",quote=F,row.names=F)
