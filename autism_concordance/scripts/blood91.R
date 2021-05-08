library(GEOquery)
library(limma)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood9.1 = getGEO(GEO= "GSE18123",GSEMatrix = TRUE,AnnotGPL = T)
blood9.1.exp = as.data.frame(exprs(blood9.1$`GSE18123-GPL6244_series_matrix.txt.gz`))
blood9.1.pdat = as.data.frame(pData(blood9.1$`GSE18123-GPL6244_series_matrix.txt.gz`))
blood9.1.info = blood9.1$`GSE18123-GPL6244_series_matrix.txt.gz`@featureData@data
blood9.1.info$`Gene symbol` = gsub("///.*","",blood9.1.info$`Gene symbol`)
rm(blood9.1)
blood9.1.exp.tmp = blood9.1.exp+abs(min(blood9.1.exp))
blood9.1.exp = log2(blood9.1.exp.tmp+1)
blood9.1.exp = normalizeQuantiles(blood9.1.exp)
blood9.1.exp = blood9.1.exp[,blood9.1.pdat$`diagnosis:ch1`=="AUTISM"|blood9.1.pdat$`diagnosis:ch1`=="CONTROL"]
blood9.1.pdat = blood9.1.pdat[blood9.1.pdat$`diagnosis:ch1`=="AUTISM"|blood9.1.pdat$`diagnosis:ch1`=="CONTROL",]
blood9.1.ed.tmp = blood9.1.pdat$`diagnosis:ch1`
blood9.1.ed.tmp = gsub("AUTISM",1,blood9.1.ed.tmp)
blood9.1.ed.tmp = gsub("CONTROL",0,blood9.1.ed.tmp)
blood9.1.ed = factor(blood9.1.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood9.1.exp.asd = blood9.1.exp[blood9.1.ed ==1]
blood9.1.exp.control = blood9.1.exp[blood9.1.ed ==0]
blood9.1.exp.control.sd = apply(blood9.1.exp.control,1,sd)
blood9.1.exp.control.mean = apply(blood9.1.exp.control,1,mean)
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
# apply the function to blood9.1 data 
blood9.1.z.df = get.z.func(blood9.1.exp.control.mean,blood9.1.exp.control.sd,blood9.1.exp.asd)
blood9.1.sig.genes = gene.name.func(blood9.1.z.df)
blood9.1.sig.gene.sym = get.gene.symbol.func(blood9.1.sig.genes,blood9.1.info)
#####################################################################################################
# separate perturbed genes by sex
blood9.1.sex = blood9.1.pdat$`gender:ch1`
blood9.1.sex = blood9.1.sex[blood9.1.pdat$`diagnosis:ch1`=="AUTISM"]
blood9.1.male.genes = blood9.1.sig.gene.sym[blood9.1.sex=="male"]
blood9.1.female.genes = blood9.1.sig.gene.sym[blood9.1.sex=="female"]
rm(blood9.1.exp.control.mean,blood9.1.exp.control.sd,blood9.1.sig.genes)

# identify differentially expressed genes via traditional t-test
blood9.1.pvals = apply(blood9.1.exp,1,function(x){t.test(x~blood9.1.ed)$p.value})
blood9.1.sig = blood9.1.pvals[blood9.1.pvals<0.05]
blood9.1.sig.names = names(blood9.1.sig)
blood9.1.sig.names = unique(blood9.1.info[blood9.1.info$ID %in% blood9.1.sig.names,]$`Gene symbol`)

