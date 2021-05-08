library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood8 = getGEO(GEO= "GSE37772",GSEMatrix = TRUE,AnnotGPL = T)
blood8.exp = as.data.frame(exprs(blood8$GSE37772_series_matrix.txt.gz))
blood8.pdat = as.data.frame(pData(blood8$GSE37772_series_matrix.txt.gz))
blood8.info = blood8$GSE37772_series_matrix.txt.gz@featureData@data
rm(blood8)
blood8.exp.tmp = blood8.exp+abs(min(blood8.exp))
blood8.exp = log2(blood8.exp.tmp+1)
blood8.exp = normalizeQuantiles(blood8.exp)
blood8.ed.tmp = blood8.pdat$`autism trait:ch1`
blood8.ed.tmp = gsub("Autism",1,blood8.ed.tmp)
blood8.ed.tmp = gsub("Control",0,blood8.ed.tmp)
blood8.ed = factor(blood8.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood8.exp.asd = blood8.exp[blood8.ed ==1]
blood8.exp.control = blood8.exp[blood8.ed ==0]
blood8.exp.control.sd = apply(blood8.exp.control,1,sd)
blood8.exp.control.mean = apply(blood8.exp.control,1,mean)
#####################################################################################################
# find which genes are perturbed in each individual
get.z.func = function(control.mean,control.sd,x){
  mean.diff = as.data.frame(apply(x,2,function(x){x - control.mean}))
  z.df = as.data.frame(apply(mean.diff,2,function(x){x/control.sd}))
  return(z.df)
}
gene.name.func = function(x){
  patient.list = list()
  for(i in 1:length(x)){
    patient = x[,i]
    sig = x[,i]>2.5 | x[,i]<(-2.5)
    sig.names=rownames(x)[sig][!is.na(rownames(x)[sig])]
    patient.list[[i]] = sig.names
  }
  return(patient.list)
}
get.gene.symbol.func=function(x,y){
  new.patient.list = list()
  for (i in 1:length(x)){
    new.info = y[x[[i]],]$`Gene symbol`
    new.info = new.info[new.info>1]
    new.patient.list[[i]] = new.info
  }
  return(new.patient.list)
}
# apply the function to blood8 data 
blood8.z.df = get.z.func(blood8.exp.control.mean,blood8.exp.control.sd,blood8.exp.asd)
blood8.sig.genes = gene.name.func(blood8.z.df)
blood8.sig.gene.sym = get.gene.symbol.func(blood8.sig.genes,blood8.info)
#####################################################################################################
# separate perturbed genes by sex
blood8.sex = blood8.pdat$`gender:ch1`
blood8.sex = blood8.sex[blood8.pdat$`autism trait:ch1` == "Autism"]
blood8.male.genes = blood8.sig.gene.sym[blood8.sex=="male"]
blood8.female.genes = blood8.sig.gene.sym[blood8.sex=="female"]
rm(blood8.exp.control.mean,blood8.exp.control.sd,blood9.sig.genes)
# identify differentially expressed genes via traditional t-test
blood8.pvals = apply(blood8.exp,1,function(x){t.test(x~blood8.ed)$p.value})
blood8.sig = blood8.pvals[blood8.pvals<0.05]
blood8.sig.names = names(blood8.sig)
blood8.sig.names = unique(blood8.info[blood8.info$ID %in% blood8.sig.names,]$`Gene symbol`)
