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
linearCox.analyze<-function(gene.code){
  
  rna<-as.numeric(rna.db[,which(colnames(rna.db)==gene.code)])
  
  #model application
  survie<-Surv(as.numeric(clinical.db$follow),clinical.db$censor)
  cox.gene<-coxph(survie ~ rna)
  
  #MEASURE TOOLS
  #aic
  AIC=AIC(cox.gene)
  
  #concordance + IC (95%)
  conc<-concordance(cox.gene)
  concordance<-conc$concordance
  concordance.IC.lower<-conc$concordance-sqrt(conc$var)*1.96
  concordance.IC.upper<-conc$concordance+sqrt(conc$var)*1.96
  
  #PLOT TOOLS
  #log(risk)
  plot(rna,log(predict(cox.gene,type="risk")),ylim=c(-1,1),ylab="log(risk)",main=paste0("log(risk) function according to RNA level of ",gene.code))
  
  #martingale residuals to RNA level
  plot(rna,residuals(cox.gene,type="martingale"),main=paste0("Martingale residuals according to RNA level of ",gene.code),ylab="martingale residuals")
  lines(smooth.spline(rna,residuals(cox.gene,type="martingale")),col='red',lwd=2.5)
  
  #Cox-Snell residuals
  resid_coxsnell <- -(residuals(cox.gene,type="martingale") - clinical.db$censor)
  fit_coxsnell <- coxph(formula = Surv(resid_coxsnell, clinical.db$censor) ~ 1,
                        ties    = c("efron","breslow","exact")[1])
  df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)
  plot(df_base_haz$time,df_base_haz$hazard,col='red',type='l',xlim=c(0,1.2),ylim=c(0,1.2),lwd=3,xlab="Cox-Snell residuals",ylab=" ",main=paste0("Cox-Snell residuals for  ",gene.code))
  abline(a=0,b=1,lty=2)

  return(list('model'=cox.gene,"AIC"=AIC,'concordance'=concordance,'concordance.IC.lower'= concordance.IC.lower,'concordance.IC.upper'=concordance.IC.upper))
}







