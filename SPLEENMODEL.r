library(Hmisc);library(glmnet);library(survivalROC);library(stats);
library(hdnom);library(boot);library(e1071);
library(reshape2);library(caret)
library(rms);library(survival);library(survcomp);library(survminer);
library(compareC);library(sva);library(MASS);
library(ggplot2);library(TSHRC);library(ggpubr)
library(survIDINRI)
library(survivalROC)
source('HLtest.r')
source('dca.r')
source('val.prob.ci.dec08.r')
source('stdca.R')
library(ggsci)
library(ggplot2)
library(gridExtra)
library(timeROC)
library(corrplot)
library(rmda)
library(ggDCA)
options (warn = -1)

R.version

data<-read.csv(file="ALL_1110222.csv",encoding="UTF-8")
data_ori <- data
head(data)
df <- data


colnames(data)



data1<- data_ori
data_ori1<-data.frame(ID=data1$ID,Center=data1$Center,group=data1$group,age=data1$age,sex=data1$sex,
                   SV=log(data1$SpleenVolume/1000),Etiology00=data1$Etiology00,Etiology0=data1$Etiology0,
                      Etiology1=data1$Etiology1,Etiology2=data1$Etiology2,
                    ALB=data1$ALB,TBIL=data1$TBIL,IBIL=data1$IBIL,DBIL=data1$DBIL,
                    ALT=data1$ALT,AST=data1$AST,SLR=data1$SLR,ALP=data1$ALP,GGT=data1$GGT,CHE=data1$CHE1,
                    BUN=data1$BUN,Cr=data1$Cr,UA=data1$UA,Hb=data1$Hb,PLT=data1$PLT,
                    INR=data1$INR,CP=data1$CP,
                    MELD=data1$MELD,ALBI=data1$ALBI,FIB=data1$FIB,ALBIFIB=data1$ALBIFIB,
                      time=data1$TIME,event=data1$Status,Status_acute=data1$Status_acute)

head(data_ori1)


###clinical and spleen volume based model building
df <- data_ori1
nje=df[which(df$Center == "NJE"),]
zj=df[which(df$Center == "ZJ"),]
zd=df[which(df$Center == "ZD"),]
wh=df[which(df$Center == "WH"),]
yz=df[which(df$Center == "YZ"),]

train_data<-rbind(zj,wh,yz)

trainx_ori<-na.omit(train_data)

trainx<-trainx_ori[,7:31]
trainy<-Surv(time=trainx_ori$time,event=trainx_ori$event)
nje_cox_glm1 <- glmnet(as.matrix(trainx),
                         trainy,
                         family = "cox",
                         type.measure = "deviance",
                         alpha=1,nfolds = nrow(trainx))
plot(nje_cox_glm1,label=T,xvar="lambda")
nje_cox_glm1 <- cv.glmnet(as.matrix(trainx),
                         trainy,
                         family = "cox",
                         type.measure = "deviance",
                         alpha=1,nfolds = nrow(trainx))
plot(nje_cox_glm1,label=T)

trainx<-trainx_ori[,6:31]
trainy<-Surv(time=trainx_ori$time,event=trainx_ori$event)
nje_cox_glm2 <- glmnet(as.matrix(trainx),
                         trainy,
                         family = "cox",
                         type.measure = "deviance",
                         alpha=1,nfolds = nrow(trainx)
                    )
plot(nje_cox_glm2,label=T,xvar="lambda")
nje_cox_glm2 <- cv.glmnet(as.matrix(trainx),
                         trainy,
                         family = "cox",
                         type.measure = "deviance",
                         alpha=1,nfolds = nrow(trainx)
                    )
plot(nje_cox_glm2,label=T)





coef1<-coef(nje_cox_glm1,s=nje_cox_glm1$lambda.1se)
active.index<-which(as.numeric(coef1)!=0)
active.coefficients<-as.numeric(coef1)[active.index]
sig_glm1<-rownames(coef1)[active.index]
sig_glm_summary1<-data.frame(coefficients=active.coefficients,name=sig_glm1)
sig_glm_summary1

coef2<-coef(nje_cox_glm2,s=nje_cox_glm2$lambda.1se)
active.index<-which(coef(nje_cox_glm2,s="lambda.1se")!=0)
active.coefficients<-as.numeric(coef2)[active.index]
sig_glm2<-rownames(coef2)[active.index]
sig_glm_summary2<-data.frame(coefficients=active.coefficients,name=sig_glm2)
sig_glm_summary2

cindex<-function(df,n,model){
    dfx<-df[,n:32]
    pre <- predict(model,as.matrix(dfx),s=c("lambda.1se"),type="link")
    dfy<- Surv(time=df$time,event=df$event)
    cdex <- rcorr.cens(pre,dfy)
    outputdf<-cbind(df,pre)
    c <- 1-cdex[1]
    upper <- cdex[3]/2*1.96+c
    lower <- c-cdex[3]/2*1.96
    CDEX<-rbind(lower,c,upper)
    print("_______")
    print(CDEX)
    
    return(pre)
}

cindex_cox <- function(inputdf,model){
    df_cox=inputdf
    surv <- Surv(time=df_cox$time,event=df_cox$event)
    pre <- predict(model,newdata =df_cox)
    outputdf<-cbind(inputdf,pre)
    cdex <- rcorr.cens(pre,surv)
    #cdex<-survConcordance(surv~pre)$concordance
    c <- 1-cdex[1]
    upper <- cdex[3]/2*1.96+c
    lower <- c-cdex[3]/2*1.96
    CDEX<-rbind(lower,c,upper)
    print("_______")
    print(CDEX)
    return(pre)
    }




train<-rbind(zj,wh,yz)
test1<-rbind(nje)
test2<-rbind(zd)

train_pred=cindex(train,6,nje_cox_glm2);train<-cbind(train,train_pred);colnames(train)[36] <-"pre_svm"
test1_pred=cindex(test1,6,nje_cox_glm2);test1<-cbind(test1,test1_pred);colnames(test1)[36] <-"pre_svm"
test2_pred=cindex(test2,6,nje_cox_glm2);test2<-cbind(test2,test2_pred);colnames(test2)[36] <-"pre_svm"

train_pred=cindex(train,7,nje_cox_glm1);train<-cbind(train,train_pred);colnames(train)[37] <-"pre_cli"
test1_pred=cindex(test1,7,nje_cox_glm1);test1<-cbind(test1,test1_pred);colnames(test1)[37] <-"pre_cli"
test2_pred=cindex(test2,7,nje_cox_glm1);test2<-cbind(test2,test2_pred);colnames(test2)[37] <-"pre_cli"

exp(log(1000))

# binding model building
## SV+CLI
sctrain <-  data.frame (time = train$time,event = train$event,pre_svm=train$pre_svm,pre_cli=train$pre_cli)
response <- Surv(time=sctrain$time,event=sctrain$event)
svclicox <- coxph(response~pre_svm+pre_cli,data=sctrain)

##sv+fib4+albi
sftrain <-  data.frame (time = train$time,event = train$event,pre_svm=train$pre_svm,ALBIFIB=train$ALBIFIB)
response <- Surv(time=sftrain$time,event=sftrain$event)
scfibcox <- coxph(response~pre_svm+ALBIFIB,data=sftrain)

##spleen volume only
svtrain <-  data.frame (time = train$time,event = train$event,spleenvolume=exp(train$SV))
response <- Surv(time=svtrain$time,event=svtrain$event)
svonecox <- coxph(response~spleenvolume,data=svtrain)



train_pred=cindex_cox(train,svclicox);train<-cbind(train,train_pred);colnames(train)[38] <-"pre_svcli"
test1_pred=cindex_cox(test1,svclicox);test1<-cbind(test1,test1_pred);colnames(test1)[38] <-"pre_svcli"
test2_pred=cindex_cox(test2,svclicox);test2<-cbind(test2,test2_pred);colnames(test2)[38] <-"pre_svcli"

train_pred=cindex_cox(train,scfibcox);train<-cbind(train,train_pred);colnames(train)[39] <-"pre_scfib"
test1_pred=cindex_cox(test1,scfibcox);test1<-cbind(test1,test1_pred);colnames(test1)[39] <-"pre_scfib"
test2_pred=cindex_cox(test2,scfibcox);test2<-cbind(test2,test2_pred);colnames(test2)[39] <-"pre_scfib"



all=rbind(train,test1,test2)
write.csv(all,"all_pre.csv",row.names = FALSE)

timeroc<- function(df,filename){
    timerocdata <- data.frame(time=df$time,status=df$event,vec1=as.vector(df$pre_svm))
    ROC<- timeROC(T=timerocdata$time, delta=timerocdata$status,
                     marker=timerocdata$vec1, cause=1,
                     weighting='cox',
                     times=c(365*3,365*5),ROC=TRUE)
    ROC_res<-ROC
    par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    plot(ROC_res,time=365*3,col="#BC3C29FF",title=FALSE,lwd=4)
    plot(ROC_res,time=365*5,col="#0072B5FF",add=TRUE,title=FALSE,lwd=4)
    
    legend('bottomright',
       c(paste0('AUC at 3 years: ',round(ROC_res$AUC[1],2)),
         paste0('AUC at 5 years: ',round(ROC_res$AUC[2],2))),
       col=c("#BC3C29FF","#0072B5FF"),lwd=5,bty = 'n')
}

find_cutoff<- function(df,filename){
    timerocdata <- data.frame(time=df$time,status=df$event,vec1=as.vector(df$pre_svm))
    predict_time1=365*3
    predict_time2=365*5
    
    auc_text=c()
    you_roc1 <- survivalROC(Stime=timerocdata$time,
                           status = timerocdata$status,
                           marker = timerocdata$vec1,
                           predict.time = predict_time1,
                           method = "KM")
    you_roc2 <- survivalROC(Stime=timerocdata$time,
                           status = timerocdata$status,
                           marker = timerocdata$vec1,
                           predict.time = predict_time2,
                           method = "KM")
    
    cutoffvalue1 <- data.frame(cutoff1=you_roc1$cut.values,
                              sensitivity_TP=you_roc1$TP,#(sensitivity)
                              specificity_FP1=you_roc1$FP,
                              cutoff2=you_roc2$cut.values,
                              sensitivity_TP1=you_roc2$TP,#(sensitivity)
                              specificity_FP2=you_roc2$FP#1-specificity#1-specificity
                             )
    write.csv(cutoffvalue1,paste0(filename,"_cutoff.csv"))
}

timeroc_acute<- function(df,filename){
    timerocdata <- data.frame(time=df$time,status=df$Status_acute,vec1=as.vector(df$pre_svm))
    ROC<- timeROC(T=timerocdata$time, delta=timerocdata$status,
                     marker=timerocdata$vec1, cause=1,
                     weighting='cox',
                     times=c(365*3,365*5),ROC=TRUE)
    ROC_res<-ROC
    #par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    #plot(ROC_res,time=365*3,col="#BC3C29FF",title=FALSE,lwd=4)
    #plot(ROC_res,time=365*5,col="#0072B5FF",add=TRUE,title=FALSE,lwd=4)
    print(
        paste0('AUC at 3 years: ',round(ROC_res$AUC[1],4),
        paste0('AUC at 5 years: ',round(ROC_res$AUC[2],4)    
    )))
    #legend('bottomright',
    #   c(paste0('AUC at 3 years: ',round(ROC_res$AUC[1],2)),
    #     paste0('AUC at 5 years: ',round(ROC_res$AUC[2],2))),
    #   col=c("#BC3C29FF","#0072B5FF"),lwd=5,bty = 'n')
}

timeroc(train,"train")
timeroc(test1,"test1")
timeroc(test2,"test2")

find_cutoff(train,"train")
find_cutoff(test1,"test1")
find_cutoff(test2,"test2")

timeroc_acute(train,"train")
timeroc_acute(test1,"test1")
timeroc_acute(test2,"test2")




calibration <- function(dcadataraw){
    dcadata <- data.frame(svm=dcadataraw$pre_svm,
                           time=dcadataraw$time,event=dcadataraw$event)
    

    data_dca <- na.omit(dcadata)
    evaluaten=floor(nrow(data_dca)/3)
    par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    traincal <- val.surv(f,S=Surv(data_dca$time,
                                  data_dca$event),
                         newdata=data_dca, u=365*3, evaluate=evaluaten)
    res = groupkm(traincal$p, Srv=Surv(data_dca$time,data_dca$event), m=evaluaten, u=365*3, pl=T, add=F,xlim=c(0,1),
                  ylim=c(0,1),errbar=T, errbar.col="#00468B",cex.axis=1,cex.lab=1,font=1,
                  xlab="Nomogram-Predicted Probability of Decompensation",lwd = 2,
                  ylab="Actual Decompensation (proportion)",cex.subtitle=F,col="#00468B")
    abline(0,1,lty=2)
    lines(res[,c('x','KM')],type= 'o',lwd = 2,col="#00468B",pch = 16)
    
    par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    traincal <- val.surv(f,S=Surv(data_dca$time,
                                  data_dca$event),
                         newdata=data_dca, u=365*5, evaluate=evaluaten)
    res = groupkm(traincal$p, Srv=Surv(data_dca$time,data_dca$event), m=evaluaten, u=365*5, pl=T, add=T,xlim=c(0,1),
                  ylim=c(0,1),errbar=T, errbar.col="#ED0000",cex.axis=1,cex.lab=1,font=1,
                  xlab="Nomogram-Predicted Probability of Decompensation",lwd = 2,
                  ylab="Actual Decompensation (proportion)",cex.subtitle=F,col="#ED0000")
    abline(0,1,lty=2)
    lines(res[,c('x','KM')],type= 'o',lwd = 2,col="#ED0000",pch = 16)
    legend('topleft',c(paste0('3-year'),paste0('5-year')),cex=1,
       col=c("#00468B","#ED0000"),lwd=5,bty = 'n')
}

dcadataraw<-train
dcadata <- data.frame(svm=dcadataraw$pre_svm,time=dcadataraw$time,event=dcadataraw$event)
f =cph(Surv(time=dcadata$time,event=dcadata$event)~svm,x = T, y = T, data  =dcadata,surv = TRUE) 


calibration(train)

calibration(rbind(test1,test2))


test_bind<-rbind(test1,test2)
df<-train
df<-test_bind
eventa<-df$event
eventa<-df$Status_acute
timerocdata1 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$pre_svm))
ROC1<- timeROC(T=timerocdata1$time, delta=timerocdata1$status,marker=timerocdata1$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)


timerocdata2 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$CP))
ROC2<- timeROC(T=timerocdata2$time, delta=timerocdata2$status,marker=timerocdata2$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata3 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$MELD))
ROC3<- timeROC(T=timerocdata3$time, delta=timerocdata3$status,marker=timerocdata3$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata4 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$FIB))
ROC4<- timeROC(T=timerocdata4$time, delta=timerocdata4$status,marker=timerocdata4$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata5 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$ALBI))
ROC5<- timeROC(T=timerocdata5$time, delta=timerocdata5$status,marker=timerocdata5$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata6 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$ALBIFIB))
ROC6<- timeROC(T=timerocdata6$time, delta=timerocdata6$status,marker=timerocdata6$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata7 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$pre_cli))
ROC7<- timeROC(T=timerocdata7$time, delta=timerocdata7$status,marker=timerocdata7$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata8 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$pre_svcli))
ROC8<- timeROC(T=timerocdata8$time, delta=timerocdata8$status,marker=timerocdata8$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)

timerocdata9 <- data.frame(time=df$time,status=eventa,vec1=as.vector(df$pre_scfib))
ROC9<- timeROC(T=timerocdata9$time, delta=timerocdata9$status,marker=timerocdata9$vec1, cause=1,weighting='cox',times=c(365*3,365*5),ROC=TRUE)


par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
plot(ROC1,time=365*3,col="#32383d",title=FALSE,lwd=3)
plot(ROC2,time=365*3,col="#be374a",add=TRUE,title=FALSE,lwd=2)
plot(ROC3,time=365*3,col="#e0e669",add=TRUE,title=FALSE,lwd=2)
plot(ROC4,time=365*3,col="#6b302f",add=TRUE,title=FALSE,lwd=2)
plot(ROC5,time=365*3,col="#88ae88",add=TRUE,title=FALSE,lwd=2)
plot(ROC6,time=365*3,col="#fe6e2d",add=TRUE,title=FALSE,lwd=2)
plot(ROC7,time=365*3,col="#32383d",add=TRUE,lty=2,title=FALSE,lwd=2)
plot(ROC8,time=365*3,col="#32383d",add=TRUE,lty=3,title=FALSE,lwd=2)
plot(ROC9,time=365*3,col="#32383d",add=TRUE,lty=4,title=FALSE,lwd=2)
#plot(ROC_res,time=365*5,col="#0072B5FF",add=TRUE,title=FALSE,lwd=4)

legend('bottomright',
   c(paste0('Spleen based model:',round(ROC1$AUC[1],2)),
     paste0('Child pugh:',round(ROC2$AUC[1],2)),
     paste0('MELD:',round(ROC3$AUC[1],2)),
     paste0('FIB-4:',round(ROC4$AUC[1],2)),
     paste0('ALBI:',round(ROC5$AUC[1],2)),
          paste0('FIB-4&ALBI: ',round(ROC6$AUC[1],2)),
     paste0('Clinical lasso model:',round(ROC7$AUC[1],2)),
     paste0('SV+CLI:',round(ROC8$AUC[1],2)),
     paste0('SV+FIB-4&ALBI:',round(ROC9$AUC[1],2))
    ),col=c("#32383d","#be374a","#e0e669","#6b302f","#88ae88","#fe6e2d","#32383d","#32383d","#32383d"),lwd=3,bty = 'n', lty=c(1,1,1,1,1,1,2,3,4),
      )
     

par(bty="o",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
plot(ROC1,time=365*5,col="#32383d",title=FALSE,lwd=3)
plot(ROC2,time=365*5,col="#be374a",add=TRUE,title=FALSE,lwd=2)
plot(ROC3,time=365*5,col="#e0e669",add=TRUE,title=FALSE,lwd=2)
plot(ROC4,time=365*5,col="#6b302f",add=TRUE,title=FALSE,lwd=2)
plot(ROC5,time=365*5,col="#88ae88",add=TRUE,title=FALSE,lwd=2)
plot(ROC6,time=365*5,col="#fe6e2d",add=TRUE,title=FALSE,lwd=2)
plot(ROC7,time=365*5,col="#32383d",add=TRUE,lty=2,title=FALSE,lwd=2)
plot(ROC8,time=365*5,col="#32383d",add=TRUE,lty=3,title=FALSE,lwd=2)
plot(ROC9,time=365*5,col="#32383d",add=TRUE,lty=4,title=FALSE,lwd=2)
#plot(ROC_res,time=365*5,col="#0072B5FF",add=TRUE,title=FALSE,lwd=4)

legend('bottomright',
   c(paste0('Spleen based model:',round(ROC1$AUC[2],2)),
     paste0('Child pugh:',round(ROC2$AUC[2],2)),
     paste0('MELD:',round(ROC3$AUC[2],2)),
     paste0('FIB-4:',round(ROC4$AUC[2],2)),
     paste0('ALBI:',round(ROC5$AUC[2],2)),
     paste0('FIB-4&ALBI: ',round(ROC6$AUC[2],2)),
     paste0('Clinical lasso model:',round(ROC7$AUC[2],2)),
     paste0('SV+CLI:',round(ROC8$AUC[2],2)),
     paste0('SV+FIB-4&ALBI:',round(ROC9$AUC[2],2))
    ),col=c("#32383d","#be374a","#e0e669","#6b302f","#88ae88","#fe6e2d","#32383d","#32383d","#32383d"),
       lwd=3,bty = 'n', lty=c(1,1,1,1,1,1,2,3,4),
      )

AUC1<-data.frame(name="pre_svm",three_year=ROC1$AUC[1],five_year=ROC1$AUC[2])
AUC2<-data.frame(name="CP",three_year=ROC2$AUC[1],five_year=ROC2$AUC[2])
AUC3<-data.frame(name="MELD",three_year=ROC3$AUC[1],five_year=ROC3$AUC[2])
AUC4<-data.frame(name="FIB-4",three_year=ROC4$AUC[1],five_year=ROC4$AUC[2])
AUC5<-data.frame(name="ALBI",three_year=ROC5$AUC[1],five_year=ROC5$AUC[2])
AUC6<-data.frame(name="ALBIFIB",three_year=ROC6$AUC[1],five_year=ROC6$AUC[2])
AUC7<-data.frame(name="pre_cli",three_year=ROC7$AUC[1],five_year=ROC7$AUC[2])
AUC8<-data.frame(name="pre_svcli",three_year=ROC8$AUC[1],five_year=ROC8$AUC[2])
AUC9<-data.frame(name="pre_scfib",three_year=ROC9$AUC[1],five_year=ROC9$AUC[2])
AUC<-rbind(AUC1,AUC2,AUC3,AUC4,AUC5,AUC6,AUC7,AUC8,AUC9)

AUC
write.csv(AUC,"AUC_atest.csv")

colnames(train)

dca_data<-train

source('stdca.R')

stdca(data=train,outcome='event',ttoutcome='time',predictors=c("pre_svm","CP","ALBI","FIB","ALBIFIB","MELD"),
      timepoint=365*3,xstop=0.55,probability=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),cmprsk=TRUE,smooth=TRUE,loess.span=0.3)


test_bind<-rbind(test1,test2)
stdca(data=test_bind,outcome='event',ttoutcome='time',predictors=c("pre_svm","CP","ALBI","FIB","ALBIFIB","MELD"),
      timepoint=365*3,xstop=0.55, probability=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),cmprsk=TRUE,smooth=TRUE,loess.span=0.3)



ccom <- function(df,filename){
    ccomapre2<-data.frame(no=1,c_sn=2,sn_lower=3,sn_upper=4,c_com=5,c_lower=6,c_upper=7,P=4,row=8)
    mydata<-df
    mydata3<-data.frame(time=mydata$time,status=mydata$event,pred=mydata$pre_svm,
                   CP=mydata$CP,MELD=mydata$MELD,
                   ALBI=mydata$ALBI,FIB=mydata$FIB,ALBIFIB=mydata$ALBIFIB,cli=mydata$pre_cli)
    i=4
    while (i<10){
        print(i)
        mydata2<-mydata3[!is.na(mydata3[,i]),]
        mydata2<-mydata2[!is.na(mydata2$pred),]
        C_index1 <- concordance.index(x=mydata2$pred, surv.time=mydata2$time,
                                    surv.event=mydata2$status,method="noether")
        C_index2 <- concordance.index(x=mydata2[,i], surv.time=mydata2$time,
                                    surv.event=mydata2$status,method="noether")
        ccomr<-cindex.comp(C_index1, C_index2)
        ccomapre1<-data.frame(no=colnames(mydata2)[i],c_sn=ccomr$cindex1,sn_lower=C_index1$lower,sn_upper=C_index1$upper,c_com=ccomr$cindex2,
                              c_lower=C_index2$lower,c_upper=C_index2$upper,P=ccomr$p.value,row=nrow(mydata2))
        ccomapre2<-rbind(ccomapre2,ccomapre1)
        i=i+1
    }
    write.csv(ccomapre2,paste0(filename,"_cindexcompare.csv"))
    }

ccom(train,"train")
ccom(test1,"test1")
ccom(test2,"test2")
ccom(rbind(test1,test2),"test_bind")

ccomacute <- function(df,filename){
    ccomapre2<-data.frame(no=1,c_sn=2,sn_lower=3,sn_upper=4,c_com=5,c_lower=6,c_upper=7,P=4,row=8)
    mydata<-df
    mydata3<-data.frame(time=mydata$time,status=mydata$Status_acute,pred=mydata$pre_svm,
                   CP=mydata$CP,MELD=mydata$MELD,
                   ALBI=mydata$ALBI,FIB=mydata$FIB,ALBIFIB=mydata$ALBIFIB,cli=mydata$pre_cli)
    i=4
    while (i<10){
        print(i)
        mydata2<-mydata3[!is.na(mydata3[,i]),]
        mydata2<-mydata2[!is.na(mydata2$pred),]
        C_index1 <- concordance.index(x=mydata2$pred, surv.time=mydata2$time,
                                    surv.event=mydata2$status,method="noether")
        C_index2 <- concordance.index(x=mydata2[,i], surv.time=mydata2$time,
                                    surv.event=mydata2$status,method="noether")
        ccomr<-cindex.comp(C_index1, C_index2)
        ccomapre1<-data.frame(no=colnames(mydata2)[i],c_sn=ccomr$cindex1,sn_lower=C_index1$lower,sn_upper=C_index1$upper,c_com=ccomr$cindex2,
                              c_lower=C_index2$lower,c_upper=C_index2$upper,P=ccomr$p.value,row=nrow(mydata2))
        ccomapre2<-rbind(ccomapre2,ccomapre1)
        i=i+1
    }
    write.csv(ccomapre2,paste0(filename,"_cindexcompare.csv"))
    }

ccomacute(train,"train_acute")
ccomacute(rbind(test1,test2),"test_bind_acute")

iniinf <- function(df,filename){
    INIDATA2<-data.frame(no=1,INI=2,UPPER=3,LOWER=4,P=5,
                     NRI=6,UPPER1=7,LOWER1=8,P1=9,row=10)
    i=4
    t0=365*5
    mydata<-df
    mydata3<-data.frame(time=mydata$time,status=mydata$event,pred=mydata$pre_svm,
                   CP=mydata$CP,MELD=mydata$MELD,
                   ALBI=mydata$ALBI,FIB=mydata$FIB,ALBIFIB=mydata$ALBIFIB,cli=mydata$pre_cli)
    #mydata2<-na.omit(mydata3)
    while (i<10){
        print(i)
        mydata2<-mydata3[!is.na(mydata3[,i]),]
        mydata2<-mydata2[!is.na(mydata2$pred),]
        indata0=mydata2[,i]
        indata1=mydata2[,3]
        covs0<-as.matrix(indata0)
        covs1<-as.matrix(indata1)
        x<-IDI.INF(mydata2[,1:2], covs0, covs1, t0, npert=1000)
        INIDATA1<-data.frame(no=colnames(mydata2)[i],INI=x$m1[1],UPPER=x$m1[2],LOWER=x$m1[3],P=x$m1[4],
                             NRI=x$m2[1],UPPER1=x$m2[2],LOWER1=x$m2[3],P1=x$m2[4],row=nrow(mydata2))
        INIDATA2<-rbind(INIDATA2,INIDATA1)
        i=i+1
    }
    write.csv(INIDATA2,paste0(filename,"_iniinf1000.csv"))
}

iniinf(train,"train")
iniinf(rbind(test1,test2),"test_bind")

#data_km<-read.csv(file="all_pre.csv",encoding="UTF-8")
#head(data_km)
#nje_km=data_km[which(data_km$Center == "NJE"),]
#zj_km=data_km[which(data_km$Center == "ZJ"),]
#zd_km=data_km[which(data_km$Center == "ZD"),]
#wh_km=data_km[which(data_km$Center == "WH"),]
#yz_km=data_km[which(data_km$Center == "YZ"),]
#train_km<-rbind(nje_km,wh_km,yz_km)
#test1_km<-zj_km
#test2_km<-zd_km

kmplot<-function(plot,n){
    df<-data.frame(pre=plot$pre_svm,time=plot$time,event=plot$event)
    df<- na.omit(df)
    veclie<- df$pre
    cutoff<- 4.4
    if(cutoff>=1){
      veclie[which(veclie < cutoff)] <- 0
      veclie[which(veclie >= cutoff)] <- 1
    }
    veclie <-  as.factor(veclie)
    #hrmatrix <- hazard.ratio(x=veclie,surv.time = df$time,surv.event = df$event)
    #hrdata1 <- data.frame(name= filename, hr = hrmatrix$hazard.ratio, low = hrmatrix$lower, up = hrmatrix$upper)
    #print(hrdata1)
    veclie1 <- veclie
    survy <- Surv(time=df$time,event = df$event)
    survkm <- survy
    
    kmdata <- data.frame(surv = survkm,vect1 = veclie1)
    kmmodel <- survfit(surv~vect1,data=kmdata)
    par(bty="l",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    ggsurv <- ggsurvplot(
        fit=kmmodel,data=kmdata,  size = 1.5, risk.table = T,     
        pval = TRUE,  conf.int = F, xlab = "Time (days)", fun = "event", xlim=c(0,n),palette = c( "lancet"),legend.labs = 
        c("Low risk","High risk"),break.time.by =365,censor=T,
        surv.plot.height=0.75,risk.table.height=0.25,ggtheme=theme_classic(),tables.theme=theme_void()
    )
    #ggsurv<- ggsurv$plot
    #ggsurv <- ggpar(
    #    ggsurv, font.title = c(12, "bold", "black"),
    #    font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
    #    font.xtickslab = c(14, "plain", "black"),  legend = "top",
    #    font.ytickslab = c(14, "plain", "black"))
    #survp <- ggsurv
    print(ggsurv)
    #summary(kmmodel)
}


kmplot(train,365*6)
kmplot(test1,365*6)
kmplot(test2,365*6)

kmplot2<-function(plot,n){
    df<-data.frame(pre=plot$pre_svm,time=plot$time,event=plot$Status_acute)
    df<- na.omit(df)
    veclie<- df$pre
    cutoff<- 4.4
    if(cutoff>=1){
      veclie[which(veclie < cutoff)] <- 0
      veclie[which(veclie >= cutoff)] <- 1
    }
    veclie <-  as.factor(veclie)
    #hrmatrix <- hazard.ratio(x=veclie,surv.time = df$time,surv.event = df$event)
    #hrdata1 <- data.frame(name= filename, hr = hrmatrix$hazard.ratio, low = hrmatrix$lower, up = hrmatrix$upper)
    #print(hrdata1)
    veclie1 <- veclie
    survy <- Surv(time=df$time,event = df$event)
    survkm <- survy
    
    kmdata <- data.frame(surv = survkm,vect1 = veclie1)
    kmmodel <- survfit(surv~vect1,data=kmdata)
    par(bty="l",pty="s", font =2,font.axis=2,font.lab=2,mfrow=c(1,1),lty=1,col=1,lwd=2)
    ggsurv <- ggsurvplot(
        fit=kmmodel,data=kmdata,  size = 1.5, risk.table = T,     
        pval = TRUE,  conf.int = F, xlab = "Time (days)", fun = "event", xlim=c(0,n),palette = c( "lancet"),legend.labs = 
        c("Low risk","High risk"),break.time.by =365,censor=T,
        surv.plot.height=0.75,risk.table.height=0.25,ggtheme=theme_classic(),tables.theme=theme_void()
    )
    #ggsurv<- ggsurv$plot
    #ggsurv <- ggpar(
    #    ggsurv, font.title = c(12, "bold", "black"),
    #    font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
    #    font.xtickslab = c(14, "plain", "black"),  legend = "top",
    #    font.ytickslab = c(14, "plain", "black"))
    #survp <- ggsurv
    print(ggsurv)
    #summary(kmmodel)
}

kmplot2(train,365*6)
kmplot2(test1,365*6)
kmplot2(test2,365*6)



datagrps<-read.csv(file="ALL_1110222Table1.csv",encoding="UTF-8")
library(CBCgrps)



head(datagrps)
colnames(datagrps)

attach(datagrps)
baseline=data.frame(group,TIME,Status,Etiology0,Etiology00,Etiology1,Etiology2,Etiology3,age, sex,ALB,TBIL,weight_kg,height_Cm,
                    IBIL,DBIL,ALT,AST,ALP,GGT,CHE,BUN,Cr,UA,Hb,PLT,INR,
                    SpleenVolume,CP,MELD,ALBI,FIB,ALBIFIB,BLEED,HE,ASCITE,HCC,Status_acute
)
varlist=c("TIME","Etiology0","Etiology00","Etiology1","Etiology2","Etiology3","age", "sex","ALB","TBIL","weight_kg","height_Cm",
          "IBIL","DBIL","ALT","AST","ALP","GGT","CHE","BUN","Cr","UA",
          "Hb","PLT","INR","SpleenVolume","CP","MELD","ALBI",	"FIB","ALBIFIB","BLEED","HE","ASCITE","HCC","Status","Status_acute"
)
tabVarlist=multigrps(baseline, gvar = "group", varlist = varlist,norm.rd = 1)
print(tabVarlist$Table, quote = T)
write.csv(tabVarlist,"TABLE1.csv")




