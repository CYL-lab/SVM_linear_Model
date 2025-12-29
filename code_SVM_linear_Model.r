

 ########################################
 #   construct the svm-forest machine   #
 ########################################


  rm(list=ls())


  # Determine the features by setting biomarker.name
  biomarker.name=c("Age", "CK18", "MGB", "WBC","Platelet")

  #############################################################################################



  list.of.packages = c("e1071", "missForest", "pROC", "ROCR")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

  library(missForest)    
  library(e1071)      
  library(pROC)
  library(ROCR)


  mrep=1000 
  bestmod_list=list()
  svm.linear=matrix(NA, mrep, 5)
  colnames(svm.linear)=c("ii", "auc", "CI.L", "CI.U", "error.rate")

  read.data=read.csv("read.data.csv",fileEncoding='big5')
  ID=c(1:398)
  read.data=cbind(read.data,ID)
  X.p=length(biomarker.name)

  input.data=read.data[, c("Name", biomarker.name, "Sample.Type","ID")]
  y=as.factor((read.data[, "Sample.Type"]=="Cancer")*1)                   # 1: Cancer; 0: Benign/Healthy


  original.data=cbind(input.data, y)                                    
  healthy=which(original.data[,'Sample.Type']=="Benign" |original.data[,'Sample.Type']=='Healthy')
  Cancer=which(original.data[,'Sample.Type']=="Cancer" )

  a=table(original.data[,'Sample.Type'])
  Healthy_proportion=(a[1]+a[3])/sum(a)
  Cancer_proportion=(a[2])/sum(a)

  set.seed(100)
  seed.mrep=100
  
  mm=48                                                                   # size of validation data     

  test_size_1=round(mm*Healthy_proportion,0)
  test_size_2=mm-test_size_1
  Healthy_test=sample(healthy, test_size_1)
  Cancer_test=sample(Cancer, test_size_2)
  test.data=rbind(original.data[Healthy_test,], original.data[Cancer_test,])
  write.table(test.data, file="test_data.csv", sep=",", row.names=F, na = "NA", fileEncoding='big5')
  test.ID=test.data[,'ID']
  svm_data=original.data[-test.ID,]

  N.total=dim(original.data)[1]
  n=N.total-mm
  train.size=round(n*3/4, 0)      
  imp.train.X_1000=array(NA, dim=c(mrep, train.size, X.p)) 
  dimnames(imp.train.X_1000)=list(rep('SVM.linear', mrep), c(1:train.size), biomarker.name)                                


  ii=1


  for (ii in 1: mrep) {

    cat("ii=", ii, "\n")
    seed.mrep=seed.mrep+1
    set.seed(seed.mrep)
  
    train_index = sample(nrow(svm_data),  train.size) 
  
    train.data = svm_data[train_index, ]
    validation.data = svm_data[-train_index, ]
    train.X=train.data[, biomarker.name]
    train.y=as.factor(train.data[, "y"])
    validation.X=validation.data[, biomarker.name]
    validation.y=as.factor(validation.data[, "y"])
  
    ##################################################################################################################
    #  Imput missing data. When there is no missing data, the imputation approach will not change the input dataset. #
    ##################################################################################################################
     # 1) impute train 
    imp.train.X = missForest(train.X)$ximp  
    imp.train.data=cbind(imp.train.X, train.y)
    imp.train.X_1000[ii,,]=as.matrix(imp.train.X)
    
    
    train.validation.X = rbind(validation.X, imp.train.X)
    imp.validation.X = missForest(train.validation.X)$ximp[1:nrow(validation.X), ] 
    imp.validation.data=cbind(imp.validation.X, validation.y)

  

    #######################################
    #  Build model SVM with linear kernel #
    #######################################
  
    tune.out=tune.svm(train.y~., data = imp.train.data,kernel='linear',
                    range=list(cost=c(0.01, 0.1 ,1 ,10 ,100 ,1000)), 
                    tunecontrol = tune.control(sampling = "cross", cross=10), 
                    probability=TRUE)
    tune.out.linear=tune.out
    bestmod=tune.out.linear$best.model     
    bestmod_list[[ii]]=tune.out.linear$best.model 
  
  
    ###########################################
    # Calculate ROC & AUC for validation data #
    ###########################################
    pred.table=table(true=validation.data[,"y"], pred=predict(tune.out.linear$best.model, imp.validation.data[, 1:X.p], probability=TRUE, decision.values =TRUE))  
  
    pred.validationdata=predict(tune.out.linear$best.model, imp.validation.data, probability=TRUE, decision.values =TRUE, fitted=TRUE)
    error.rate.linear=(sum(pred.table)-sum(diag(pred.table)))/dim(imp.validation.data)[1]
    cat("svm.error.rate based on linear kernel=", error.rate.linear, "\n")
    print(pred.table)
   
  
    pred_validation = attributes(predict(tune.out.linear$best.model, imp.validation.data, probability=TRUE, decision.values =TRUE, fitted=TRUE))$prob
  
    roc.results=roc(validation.y, pred_validation[, "1"], ci=TRUE)   
    roc.out=cbind(roc.results$sensitivities, roc.results$specificities, roc.results$thresholds)
    colnames(roc.out)=c("sensitivity", "specificity", "threshold")
    plot(1-roc.out[,"specificity"], roc.out[,"sensitivity"], type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), 
       main="Validation data by SVM with linear kernel \n with imputation")
    lines(1-roc.out[,"specificity"], roc.out[,"sensitivity"], lwd=2)
    abline(0, 1, lty=2)
    text(0.7, 0.3, paste("AUC=", round(roc.results$auc, 2)))
    text(0.7, 0.2, paste("95 % CI =", round(roc.results$ci[1],2), "-", round(roc.results$ci[3],2)))
  
    svm.linear[ii, ]=c(ii, roc.results$auc, as.numeric(roc.results$ci)[1], as.numeric(roc.results$ci)[3], error.rate.linear)
  

  
  } 


  cat("linear:", "\n"); print(round(apply(svm.linear, 2, mean), 2)[-1])


  cat("linear(sd):", "\n"); print(round(apply(svm.linear, 2, sd), 2)[-1])

  write.table(svm.linear, "svm.linear(model1).txt", quote=F)
    

  saveRDS(bestmod_list, file = "svm_bestmodel.rds", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  saveRDS(imp.train.X_1000, file="imp.train.X_1000.rds", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)



