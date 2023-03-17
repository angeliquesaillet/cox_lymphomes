###############################################################################

# 4. Automated procedure to be applied on each gene (only first 4 for the moment)

#Author : Angélique Saillet 
#mail : angelique.saillet@etu.univ-grenoble-alpes.fr
############################################################################## 

source("../R/mutations_significance.R") #loading the previous file to recover the datasets

#function that applies the automated survival analye according to mutations effect on the RNA and linear cox performance(adequate or not)
automated.survival.analyze<-function(gene.code,graph=TRUE){
  
  if(signif.mutations.table[which(signif.mutations.table$gene==gene.code),'test']==1){
    return(linearCox.analyze(gene.code,mutations=F,graph=graph))#only with the non-mutated observations
  }else{
    return(linearCox.analyze(gene.code,mutations="ALL",graph=graph))#with all the observations, regardless the mutation info.
  }
  
  #decide if the linear cox is satisfactory with the help of previous output
  
  ## if not we apply the non-linear cox model 
  
  
}