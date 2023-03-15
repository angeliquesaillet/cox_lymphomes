###############################################################################

# 2. Recovery of the significativity of the mutation effect for each gene

#Author : Angélique Saillet 
#mail : angelique.saillet@etu.univ-grenoble-alpes.fr
############################################################################## 

source("../R/linear_cox.R") #loading the previous file to recover the datasets

#------------------------------------------------------------------------------
#Identification des différences des ensembles de gènes entre RNA et mutations

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

#3 gènes de mutations ne sont pas contenus dans le dataset RNA : on les retire donc de mutations

mutations.db<-mutations.db[,-idc.diff]
#------------------------------------------------------------------------------

#Function that calculate an significativity indicator for each genes according to the rna and mutations informations
significativity.indicator.table<-function(){

pv.test<-c()  
pv.cox.mutation<-c()
pv.cox.interaction<-c()

test<-rep(0,length(colnames(mutations.db)[-1]))
cox.model<-rep(0,length(colnames(mutations.db)[-1]))
i=1
for (gene.code in colnames(mutations.db)[-1]){
# On récupere les observations rna en fonction de leur statut de mutation  et les mutations
rna<-as.numeric(rna.db[,which(colnames(rna.db)==gene.code)])
idc.m<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==1)
rna.m<-rna[idc.m]
idc.nm<-which(mutations.db[,which(colnames(mutations.db)==gene.code)]==0)
rna.nonm<-rna[idc.nm]
mutations<-unlist(mutations.db[,which(colnames(mutations.db)==gene.code)])


#1. tester si l’expression des gènes est similaire dans les cas mutés/non mutés (tests de comparaison de moyennes)
### test de normalité préliminaire, H0 = hypothèse de normalité
shap.m.pv<-shapiro.test(rna.m)$p.value
shap.nm.pv<-shapiro.test(rna.nonm)$p.value

if((shap.m.pv>0.05)&(shap.nm.pv>0.05)){
  ### test de variance 
  var.test.pv<-var.test(rna.m,rna.nonm)$p.value  # H0 : variance égales & H1 : variance différente
  if(var.test.pv>0.05){
    ### test de student 
    pv.test<-c(pv.test,t.test(rna.m,rna.nonm,var.equal = T)$p.value) #H0 : µ1 = µ2 – H1 : µ1 ≠ µ2 (bilatérale) 
  }else{
    pv.test<-c(pv.test,t.test(rna.m,rna.nonm,var.equal = F)$p.value) #H0 : µ1 = µ2 – H1 : µ1 ≠ µ2 (bilatérale) 
  }
}else{
  pv.test<-c(pv.test,wilcox.test(rna.m,rna.nonm)$p.value)
}
#si pv < 0.05 alors les moyennes ne sont pas égales et donc effet de la mutation 


#2. Intégration de la mutation (et son interaction) dans le modèle de COX en tant que variable 
#model application
survie<-Surv(as.numeric(clinical.db$follow),clinical.db$censor)
cox.gene<-coxph(survie ~ rna*mutations)
pv.cox.mutation<-c(pv.cox.mutation,summary(cox.gene)$coefficients[2,5])
pv.cox.interaction<-c(pv.cox.interaction,summary(cox.gene)$coefficients[3,5])

#Potentiel ajout d'un ajustement des p.valeurs 

if(pv.test[i]<0.05){
  test[i]<-1
}
if((pv.cox.mutation[i]<0.05)||(pv.cox.interaction[i]<0.05)){
  cox.model[i]<-1
}
i=i+1
}
# En sortie : On doit récupérer une liste contenant un indicateur binaire de significativité pour chaque gène (dont on a l'info de mutation)

#Dataset building
results.table=data.frame(cbind(colnames(mutations.db)[-1],test,cox.model))
colnames(results.table)[1]<-'gene'

return(results.table)
}



