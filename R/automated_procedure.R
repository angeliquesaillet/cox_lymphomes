###############################################################################

# 4. Automated procedure to be applied on each gene (only first 4 for the moment)

#Author : Angélique Saillet 
#mail : angelique.saillet@etu.univ-grenoble-alpes.fr
############################################################################## 

source("../R/mutations_significance.R") #loading the previous file to recover the datasets

#function that applies the automated survival analye according to mutations effect on the RNA, [add if more]
automated.survival.analyze<-function(gene.code){
  
  if(signif.mutations.table[which(signif.mutations.table$gene==gene.code),'test']==1){
    linearCox.analyze(gene.code,mutations=F,graph=TRUE)#only with the non-mutated observations
  }else{
    linearCox.analyze(gene.code,mutations="ALL",graph=TRUE)#with all the observations, regardless the mutation info.
  }
  
}