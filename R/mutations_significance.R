###############################################################################

# 2. Recovery of the significance of the mutation effect for each gene

#Author : Angélique Saillet 
#mail : angelique.saillet@etu.univ-grenoble-alpes.fr
############################################################################## 

source("../R/linear_cox.R") #loading the previous file to recover the datasets

#------------------------------------------------------------------------------
#Identification of gene set differences between RNA and mutations

genes.mutations<-colnames(mutations.db)[-1]
rna.gene<-colnames(rna.db)[-1]
genes.diff<-c()
idc.diff<-c()

for(i in 1:length(genes.mutations)){
  if(sum(genes.mutations[i]==rna.gene)==0){
    genes.diff<-c(genes.diff,genes.mutations[i])
    idc.diff<-c(idc.diff,i)
  }
}
idc.diff<-idc.diff+1

#3 genes from the mutation dataset are not contained in the RNA dataset: they are therefore removed from the mutation dataset
mutations.db<-mutations.db[,-idc.diff]

#------------------------------------------------------------------------------

#Function that calculate a significance indicator for each genes according to the rna and mutations information
significance.indicator.table<-function(){
  #the function returns a list containing a binary indicator of significance 
  #for each gene (for which we have the mutation info i.e. 147 genes)
  
  pv.test<-c()  
  pv.cox.mutation<-c()
  pv.cox.interaction<-c()
  
  test<-rep(0,length(colnames(mutations.db)[-1]))
  cox.model<-rep(0,length(colnames(mutations.db)[-1]))
  i=1
  for (gene.code in colnames(mutations.db)[-1]){
    # We divide the rna observations according to their mutation status 
    rna<-as.numeric(rna.db[,which(colnames(rna.db)==gene.code)])
    idc.m<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==1)
    rna.m<-rna[idc.m]
    idc.nm<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==0)
    rna.nonm<-rna[idc.nm]
    mutations<-unlist(mutations.db[,which(colnames(mutations.db)==gene.code)])
    
    
    #1. testing whether the expression of the genes is similar in the mutated/non-mutated cases (student,wilcoxon & welch test)
    ### Shapiro test H0 = noramlity
    shap.m.pv<-shapiro.test(rna.m)$p.value
    shap.nm.pv<-shapiro.test(rna.nonm)$p.value
    
    if((shap.m.pv>0.05)&(shap.nm.pv>0.05)){
      ### variance test
      var.test.pv<-var.test(rna.m,rna.nonm)$p.value  # H0 : variance égales & H1 : variance différente
      if(var.test.pv>0.05){
        ### student test
        pv.test<-c(pv.test,t.test(rna.m,rna.nonm,var.equal = T)$p.value) #H0 : µ1 = µ2 – H1 : µ1 ≠ µ2 (bilatérale) 
      }else{
        ### welch test
        pv.test<-c(pv.test,t.test(rna.m,rna.nonm,var.equal = F)$p.value) #H0 : µ1 = µ2 – H1 : µ1 ≠ µ2 (bilatérale) 
      }
    }else{
      ### wilcoxon test
      pv.test<-c(pv.test,wilcox.test(rna.m,rna.nonm)$p.value)
    }
    # INTERPRETATION : if pv < 0.05 then the means are not equal and therefore effect of the mutation 
    
    #2. Integration of the mutation (and its interaction) in the COX model as a variable 
    #model application
    survie<-Surv(as.numeric(clinical.db$follow),clinical.db$censor)
    cox.gene<-coxph(survie ~ rna*mutations)
    pv.cox.mutation<-c(pv.cox.mutation,summary(cox.gene)$coefficients[2,5])
    pv.cox.interaction<-c(pv.cox.interaction,summary(cox.gene)$coefficients[3,5])
    
    ### HERE : Potential addition of a p.values adjustment ###
    pv.test.adjusted<-p.adjust(pv.test,"holm")
    pv.cox.mutation.adjusted<-p.adjust(pv.cox.mutation,"holm")
    pv.cox.interaction.adjusted<-p.adjust(pv.cox.interaction,"holm")
    
    
    
    #Determination of the significance indicator
    if(pv.test[i]<0.05){
      test[i]<-1
    }
    if((pv.cox.mutation[i]<0.05)||(pv.cox.interaction[i]<0.05)){
      cox.model[i]<-1
    }
    i=i+1
  }
  
  #Results dataset building
  results.table=data.frame(cbind(colnames(mutations.db)[-1],test,cox.model,pv.test,pv.test.adjusted,pv.cox.mutation,pv.cox.mutation.adjusted,pv.cox.interaction,pv.cox.interaction.adjusted))
  colnames(results.table)[1]<-'gene'
  
  return(results.table)
}

signif.mutations.table<-significance.indicator.table()


