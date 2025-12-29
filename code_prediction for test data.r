


  ##################################
  #    svm_forest for test data    #
  ##################################



  rm(list=ls())

  list.of.packages = c("e1071", "missForest", "pROC", "ROCR")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

  library(missForest)    
  library(e1071)      
  library(pROC)
  library(ROCR)




  bestmod_list=readRDS(file = "svm_bestmodel.rds")
  imp.train.X_1000=readRDS(file = "imp.train.X_1000.rds")
  test.data=read.csv("test_data.csv",fileEncoding='big5')

  mrep=length(bestmod_list)
  mm=nrow(test.data)
  biomarker.name=c("Age", "CK18", "MGB", "WBC","Platelet")

  pred_test.out=array(NA, dim=c(mrep, mm, 3))    
  voting.out=matrix(NA, mm, 2)
  probability=matrix(NA,mm, 1)
  colnames(voting.out)=c("Voting.y", "True.y")

  forest_data.X=test.data[, biomarker.name]
  forest_data.y=test.data[, "y"]

  for (i in 1:mrep) {
      forest.train.X = rbind(forest_data.X, imp.train.X_1000[i,,])
      imp.forest.X = missForest(forest.train.X)$ximp[1:nrow(forest_data.X), ] 
      imp.forest_data=cbind(imp.forest.X, forest_data.y)

      out= attributes(predict(bestmod_list[[i]], imp.forest_data, probability=TRUE, decision.values =TRUE, fitted=TRUE))$prob
      pred_test.out[i, ,]=cbind(rep(i, mm),  out[,'1'],  (out[,'1']>.5)*1)   
  }

  for (k in 1:mm) {
    voting.out[k,]=c((apply(pred_test.out[,k,], 2, mean)[3]>=0.5)*1 , as.numeric(imp.forest_data[k, "forest_data.y"])-1)
    probability[k,]=apply(pred_test.out[,k,],2,mean)[2]  
  }


  table.pred=table(voting.out[,1], voting.out[,2])
  TP=table.pred[2,2]
  FP=table.pred[2,1]
  FN=table.pred[1,2]
  precision=TP/(TP+FP)
  recall=TP/(TP+FN)  
  specificity=mean(1-voting.out[voting.out[, "True.y"]==0, "Voting.y"])
  F1_score=(2*precision*recall)/(precision+recall)


  cat("Votint results=", "\n")
  print(table.pred)

  cat("precision=", precision, "\n")
  cat("Sensitivity=", recall, "\n")
  cat("Specificity=", specificity, "\n")
  cat("F1_score=", F1_score, "\n")

  cat("accuracy rate=", sum(diag(table(voting.out[,1], voting.out[,2])))/mm, "\n")
  cat("Forest of svm predicted probability=", "\n"); print(apply(pred_test.out[,,2], 2, mean))


  new.data=cbind(test.data[, "y"], voting.out,probability)

  colnames(new.data)=c('Stage',"Prediction",'True','Probability')



  ########################
  #  Draw the ROC curve  #
  ########################
  roc.results=roc(new.data[, "True"], new.data[, "Probability"], ci=TRUE)   
  roc.out=cbind(roc.results$sensitivities, roc.results$specificities, roc.results$thresholds)
  colnames(roc.out)=c("sensitivity", "specificity", "threshold")
  plot(1-roc.out[,"specificity"], roc.out[,"sensitivity"], type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), 
       main="Validation data by SVM with linear kernel \n with imputation")
  lines(1-roc.out[,"specificity"], roc.out[,"sensitivity"], lwd=2)
  abline(0, 1, lty=2)
  text(0.7, 0.3, paste("AUC=", round(roc.results$auc, 2)))
  text(0.7, 0.2, paste("95 % CI =", round(roc.results$ci[1],2), "-", round(roc.results$ci[3],2)))
  

 sink("sessionInfo.txt")
 sessionInfo()
 sink()

