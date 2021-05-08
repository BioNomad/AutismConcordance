library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood12 = getGEO(GEO= "GSE25507",GSEMatrix = TRUE,AnnotGPL = T)
blood12.exp = as.data.frame(exprs(blood12$GSE25507_series_matrix.txt.gz))
blood12.pdat = as.data.frame(pData(blood12$GSE25507_series_matrix.txt.gz))
blood12.info = blood12$GSE25507_series_matrix.txt.gz@featureData@data
blood12.info$`Gene symbol` = gsub("///.*","",blood12.info$`Gene symbol`)
rm(blood12)
blood12.exp.tmp = blood12.exp+abs(min(blood12.exp))
blood12.exp = log2(blood12.exp.tmp+1)
blood12.exp = normalizeQuantiles(blood12.exp)
blood12.ed.tmp = blood12.pdat$`group:ch1`
blood12.ed.tmp = gsub("autism",1,blood12.ed.tmp)
blood12.ed.tmp = gsub("control",0,blood12.ed.tmp)
blood12.ed = factor(blood12.ed.tmp,levels = c(0,1))
######################################################################################################
# identify normal expression range based on control group 
blood12.exp.asd = blood12.exp[blood12.ed ==1]
blood12.exp.control = blood12.exp[blood12.ed ==0]
blood12.exp.control.sd = apply(blood12.exp.control,1,sd)
blood12.exp.control.mean = apply(blood12.exp.control,1,mean)
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
blood12.z.df = get.z.func(blood12.exp.control.mean,blood12.exp.control.sd,blood12.exp.asd)
blood12.sig.genes = gene.name.func(blood12.z.df)
blood12.sig.gene.sym = get.gene.symbol.func(blood12.sig.genes,blood12.info)
#####################################################################################################
rm(blood12.exp.control.mean,blood12.exp.control.sd,blood12.sig.genes)
# identify differentially expressed genes via traditional t-test
blood12.pvals = apply(blood12.exp,1,function(x){t.test(x~blood12.ed)$p.value})
blood12.sig = blood12.pvals[blood12.pvals<0.05]
blood12.sig.names = names(blood12.sig)
blood12.sig.names = unique(blood12.info[blood12.info$ID %in% blood12.sig.names,]$`Gene symbol`)
