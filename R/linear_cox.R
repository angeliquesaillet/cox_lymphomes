###############################################################################

# 1. data pre-treatment & Cox linear function

#Author : Angélique Saillet 
#mail : angelique.saillet@etu.univ-grenoble-alpes.fr
############################################################################## 


# function to treate NAs & gene matching from the the 2 data bases (if 3 -> 2 calls)
#return the db without NAs who contains same id set in the right order
db.pretreatment<-function(db.clinical,db.rna){
  #Numeric transformation if necessary
  db.clinical$follow<-as.numeric(db.clinical$follow)
  
  #id set & order
  idc.C<-c()
  idc.G<-c()
  for(i in 1:length(db.clinical$id)){
    if(length(which(db.rna$id==db.clinical$id[i]))>0){
      idc.C<-c(idc.C,i)
      idc.G<-c(idc.G,which(db.rna$id==db.clinical$id[i]))
    }
  }
  
  db.clinical<-db.clinical[idc.C,]
  db.rna<-db.rna[idc.G,]
  
  #NAs
  if(length(which(is.na(db.clinical$follow)))+length(which(is.na(db.clinical$censor)))>0){
  idc.NA<-unique(c(which(is.na(db.clinical$follow)),which(is.na(db.clinical$censor))))
  db.clinical<-db.clinical[-idc.NA,]
  db.rna<-db.rna[-idc.NA,]
  }
  return(list("clinical.db"=db.clinical,'rna.db'=db.rna))
}

#--------------------------------------------------------------
#Load original datasets 
library(readxl)
clinical.original <- read_excel("../Data/Duke_Follow_censor_ABC_GCB_Hans.xlsx")
rna.original<-read.csv2("../Data/original data/Duke-1001DLBCL_log2_rpm.csv",sep=",")
colnames(rna.original)[1]="id"
mutation.original <- read_excel("../Data/Duke_Mutations.xlsx")


#Transform datasets by using the function 
clinical.db<-db.pretreatment(clinical.original,rna.original)$clinical.db
rna.db<-db.pretreatment(clinical.original,rna.original)$rna.db
mutations.db<-db.pretreatment(clinical.db,mutation.original)$rna.db
#--------------------------------------------------------------

#function that apply a Cox linear-gene model 
#return different output as the model, R2, AIC & different plots (ROC, martingale,Cox Snell & Score residues,Concordance with IC)
library(survival)
library(survivalROC)
library(CoxR2)

linearCox.analyze<-function(gene.code,mutations="ALL",graph=TRUE){
  # the parameter gene.code defines the specific gene on which the model is applied
  
  #the mutations parameter allows to decide on which data set the model is applied according to the mutation status of the observations i.e. 
  ## if mutations = "ALL" then the model is applied on all observations independently of their mutation status (by default)
  ## if mutations = TRUE then the model is applied only on the observations for which the gene is mutated. 
  ## if mutations = FALSE then the model is applied only on the observations for which the gene is not mutated. 
  
  #The graph parameter allows to decide if we want to show the graphics results or not, TRUE by default. 
  
  #Observations selection according to the mutations parameter
  rna<-as.numeric(rna.db[,which(colnames(rna.db)==gene.code)])
  clinical<-clinical.db
  
  if(mutations==T){
    idc.m<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==1)
    rna<-rna[idc.m]
    clinical<-clinical[idc.m,]
  }
  
  if(mutations==F){
    idc.nm<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==0)
    rna<-rna[idc.nm]
    clinical<-clinical[idc.nm,]
  }
  
  
  #model application
  survie<-Surv(as.numeric(clinical$follow),clinical$censor)
  cox.gene<-coxph(survie ~ rna)
  
  #MEASURE TOOLS
  #r2
  R.squared<-coxr2(cox.gene)$rsq 
  
  #aic
  AIC=AIC(cox.gene)
  
  #concordance + IC (95%)
  conc<-concordance(cox.gene)
  concordance<-conc$concordance
  concordance.IC.lower<-conc$concordance-sqrt(conc$var)*1.96
  concordance.IC.upper<-conc$concordance+sqrt(conc$var)*1.96
  
  #P.value d'un test de rapport de vraisemblance 
  LRatio.pvalue<-summary(cox.gene)$logtest[3]
  
  #P.value du test de Log rank 
  LRank.pvalue<-summary(cox.gene)$sctest[3]
  
  #the tools used to make the titles indicative 
  if(mutations=="ALL"){
    mutations.title="all observations"
  }
  if(mutations==T){
    mutations.title="only with the mutated observations"
  }
  if(mutations==F){
    mutations.title="only with the non-mutated observations"
  }
  
  if(graph==T){
    #PLOT TOOLS
    #log(risk)
    plot(rna,log(predict(cox.gene,type="risk")),ylim=c(-1,1),ylab="log(risk)",main=paste0("log(risk) function according to RNA level of ",gene.code,"\n *",mutations.title,"*"))
    
    #martingale residuals to RNA level
    plot(rna,residuals(cox.gene,type="martingale"),main=paste0("Martingale residuals according to RNA level of ",gene.code),ylab="martingale residuals")
    lines(smooth.spline(rna,residuals(cox.gene,type="martingale")),col='red',lwd=2.5)
    
    #Cox-Snell residuals
    resid_coxsnell <- -(residuals(cox.gene,type="martingale") - clinical$censor)
    fit_coxsnell <- coxph(formula = Surv(resid_coxsnell, clinical$censor) ~ 1,
                          ties    = c("efron","breslow","exact")[1])
    df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)
    plot(df_base_haz$time,df_base_haz$hazard,col='red',type='l',xlim=c(0,max(df_base_haz$time)+0.2),ylim=c(0,max(df_base_haz$hazard)+0.2),lwd=3,xlab="Cox-Snell residuals",ylab=" ",main=paste0("Cox-Snell residuals for  ",gene.code))
    abline(a=0,b=1,lty=2)
    
    #Courbes ROC 
    scr<-predict(cox.gene,type='risk')
    t<-10
    roc.1<-survivalROC(Stime=clinical$follow,status=clinical$censor,marker=scr,predict.time=t,span=0.05)
    roc.2<-survivalROC(Stime=clinical$follow,status=clinical$censor,marker=scr,predict.time=t,method="KM")
    plot(roc.1$FP,roc.1$TP,type='l',col='red',main=paste0("Courbe ROC à l'instant t = ",t),ylab="TP",xlab="FP")
    points(roc.2$FP,roc.2$TP,type='l',col='blue')
    text(c(0.8,0.8),c(0.2,0.1),labels=c(paste0("AUC (NNE) = ",round(roc.1$AUC,3)),paste0("AUC (KM) = ",round(roc.2$AUC,3))),col=c('red','blue'))
    abline(a=0,b=1)
  }
  
  return(list('model'=cox.gene,"AIC"=AIC,"concordance"=list('value'=concordance,'concordance.IC.lower'= concordance.IC.lower,'concordance.IC.upper'=concordance.IC.upper),'LRatio.pvalue'=LRatio.pvalue,'LRank.pvalue'=LRank.pvalue,'R.squared'=R.squared))
}







