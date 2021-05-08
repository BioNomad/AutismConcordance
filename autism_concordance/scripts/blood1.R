library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood1 = getGEO(GEO= "GSE26415",GSEMatrix = TRUE,AnnotGPL = T)
blood1.exp = as.data.frame(exprs(blood1$GSE26415_series_matrix.txt.gz))
blood1.pdat = as.data.frame(pData(blood1$GSE26415_series_matrix.txt.gz))
blood1.info = blood1$GSE26415_series_matrix.txt.gz@featureData@data
blood1.exp.tmp = blood1.exp+abs(min(blood1.exp))
blood1.exp = log2(blood1.exp.tmp+1)
blood1.exp = normalizeQuantiles(blood1.exp)
rm(blood1)
blood1.ed.tmp = blood1.pdat$`sample type:ch1`[blood1.pdat$`sample type:ch1`=="ASD" | 
                                                       blood1.pdat$`sample type:ch1` == "control"]
blood1.exp = blood1.exp[,blood1.pdat$`sample type:ch1`=="ASD" | 
                          blood1.pdat$`sample type:ch1` == "control"]
blood1.pdat.select = blood1.pdat[blood1.pdat$`sample type:ch1`=="ASD" | 
                                   blood1.pdat$`sample type:ch1` == "control",]
blood1.ed.tmp = gsub("ASD", 1,blood1.ed.tmp)
blood1.ed.tmp = gsub("control",0,blood1.ed.tmp)
blood1.ed = factor(blood1.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood1.exp.asd = blood1.exp[blood1.ed ==1]
blood1.exp.control = blood1.exp[blood1.ed ==0]
blood1.exp.control.sd = apply(blood1.exp.control,1,sd)
blood1.exp.control.mean = apply(blood1.exp.control,1,mean)
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
# apply the function to blood1 data 
blood1.z.df = get.z.func(blood1.exp.control.mean,blood1.exp.control.sd,blood1.exp.asd)
blood1.sig.genes = gene.name.func(blood1.z.df)
blood1.sig.gene.sym = get.gene.symbol.func(blood1.sig.genes,blood1.info)
##################################################################################################
rm(blood1.exp.control.mean,blood1.exp.control.sd,blood1.sig.genes)
#blood1.keep = apply(blood1.exp, 1, function(x) length(unique(x[!is.na(x)])) != 1)
#blood1.exp.rem = blood1.exp[blood1.keep,]
# identify differentially expressed genes via traditional t-test
blood1.pvals = apply(blood1.exp.rem,1,function(x){t.test(x~blood1.ed)$p.value})
blood1.sig = blood1.pvals[blood1.pvals<0.05]
blood1.sig.names = names(blood1.sig)
blood1.sig.names = unique(blood1.info[blood1.info$ID %in% blood1.sig.names,]$`Gene symbol`)
