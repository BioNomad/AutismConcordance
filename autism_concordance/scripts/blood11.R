library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood11 = getGEO(GEO= "GSE42133",GSEMatrix = TRUE,AnnotGPL = T)
blood11.exp = as.data.frame(exprs(blood11$GSE42133_series_matrix.txt.gz))
blood11.pdat = as.data.frame(pData(blood11$GSE42133_series_matrix.txt.gz))
blood11.info = blood11$GSE42133_series_matrix.txt.gz@featureData@data
rm(blood11)
blood11.ed.tmp = blood11.pdat$`dx (diagnosis):ch1`
blood11.ed.tmp = gsub("ASD",1,blood11.ed.tmp)
blood11.ed.tmp = gsub("Control",0,blood11.ed.tmp)
blood11.ed = factor(blood11.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood11.exp.asd = blood11.exp[blood11.ed ==1]
blood11.exp.control = blood11.exp[blood11.ed ==0]
blood11.exp.control.sd = apply(blood11.exp.control,1,sd)
blood11.exp.control.mean = apply(blood11.exp.control,1,mean)
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
# apply the function to blood11 data 
blood11.z.df = get.z.func(blood11.exp.control.mean,blood11.exp.control.sd,blood11.exp.asd)
blood11.sig.genes = gene.name.func(blood11.z.df)
blood11.sig.gene.sym = get.gene.symbol.func(blood11.sig.genes,blood11.info)
#####################################################################################################
rm(blood11.exp.control.mean,blood11.exp.control.sd,blood11.sig.genes)
# identify differentially expressed genes via traditional t-test
blood11.pvals = apply(blood11.exp,1,function(x){t.test(x~blood11.ed)$p.value})
blood11.sig = blood11.pvals[blood11.pvals<0.05]
blood11.sig.names = names(blood11.sig)
blood11.sig.names = unique(blood11.info[blood11.info$ID %in% blood11.sig.names,]$`Gene symbol`)
