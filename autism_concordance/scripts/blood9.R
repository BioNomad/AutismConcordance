library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood9 = getGEO(GEO= "GSE18123",GSEMatrix = TRUE,AnnotGPL = T)
blood9.exp = as.data.frame(exprs(blood9$`GSE18123-GPL570_series_matrix.txt.gz`))
blood9.pdat = as.data.frame(pData(blood9$`GSE18123-GPL570_series_matrix.txt.gz`))
blood9.info = blood9$`GSE18123-GPL570_series_matrix.txt.gz`@featureData@data
blood9.info$`Gene symbol`= gsub("///.*","",blood9.info$`Gene symbol`)
rm(blood9)
blood9.exp.tmp = blood9.exp+abs(min(blood9.exp))
blood9.exp = log2(blood9.exp.tmp+1)
blood9.exp = normalizeQuantiles(blood9.exp)
blood9.exp = blood9.exp[,blood9.pdat$`diagnosis:ch1`=="AUTISM"|blood9.pdat$`diagnosis:ch1`=="CONTROL"]
blood9.pdat = blood9.pdat[blood9.pdat$`diagnosis:ch1`=="AUTISM"|blood9.pdat$`diagnosis:ch1`=="CONTROL",]
blood9.ed.tmp = blood9.pdat$`diagnosis:ch1`
blood9.ed.tmp = gsub("AUTISM",1,blood9.ed.tmp)
blood9.ed.tmp = gsub("CONTROL",0,blood9.ed.tmp)
blood9.ed = factor(blood9.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood9.exp.asd = blood9.exp[blood9.ed ==1]
blood9.exp.control = blood9.exp[blood9.ed ==0]
blood9.exp.control.sd = apply(blood9.exp.control,1,sd)
blood9.exp.control.mean = apply(blood9.exp.control,1,mean)
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
# apply the function to blood9 data 
blood9.z.df = get.z.func(blood9.exp.control.mean,blood9.exp.control.sd,blood9.exp.asd)
blood9.sig.genes = gene.name.func(blood9.z.df)
blood9.sig.gene.sym = get.gene.symbol.func(blood9.sig.genes,blood9.info)
#####################################################################################################
rm(blood9.exp.control.mean,blood9.exp.control.sd,blood9.sig.genes)
# identify differentially expressed genes via traditional t-test
blood9.pvals = apply(blood9.exp,1,function(x){t.test(x~blood9.ed)$p.value})
blood9.sig = blood9.pvals[blood9.pvals<0.05]
blood9.sig.names = names(blood9.sig)
blood9.sig.names = unique(blood9.info[blood9.info$ID %in% blood9.sig.names,]$`Gene symbol`)
