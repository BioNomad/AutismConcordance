library(GEOquery)
library(limma)
library(DGCA)
######################################################################################################
# load gene expression data and normalize if needed
# load feature data for gene symbol annotations
# load patient data for separation by disease status
blood3 = getGEO(GEO= "GSE6575",GSEMatrix = TRUE,AnnotGPL = T)
blood3.exp = as.data.frame(exprs(blood3$GSE6575_series_matrix.txt.gz))
blood3.pdat = as.data.frame(pData(blood3$GSE6575_series_matrix.txt.gz))
blood3.info = blood3$GSE6575_series_matrix.txt.gz@featureData@data
blood3.info$`Gene symbol`= gsub("///.*","",blood3.info$`Gene symbol`)
#blood3.exp.tmp = blood3.exp+abs(min(blood3.exp))
#blood3.exp = log2(blood3.exp.tmp+1)
#blood3.exp = normalizeQuantiles(blood3.exp)
rm(blood3)
blood3.ed.tmp = blood3.pdat$characteristics_ch1
blood3.ed.tmp = gsub(",.*","",blood3.ed.tmp)
blood3.ed.tmp = gsub(" .*","",blood3.ed.tmp)
blood3.pdat$diagnosis = blood3.ed.tmp
blood3.exp = blood3.exp[,blood3.pdat$diagnosis=="Autism"|
                          blood3.pdat$diagnosis=="General"]
blood3.pdat.select = blood3.pdat[blood3.pdat$diagnosis=="Autism"|
                                   blood3.pdat$diagnosis=="General",]
blood3.ed.tmp = blood3.pdat$diagnosis[blood3.pdat$diagnosis=="Autism"|
                                         blood3.pdat$diagnosis=="General"]
blood3.ed.tmp = gsub("Autism",1,blood3.ed.tmp)
blood3.ed.tmp = gsub("General",0,blood3.ed.tmp)
blood3.ed = factor(blood3.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify normal expression range based on control group 
blood3.exp.asd = blood3.exp[blood3.ed ==1]
blood3.exp.control = blood3.exp[blood3.ed ==0]
blood3.exp.control.sd = apply(blood3.exp.control,1,sd)
blood3.exp.control.mean = apply(blood3.exp.control,1,mean)
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
# apply the function to blood3 data 
blood3.z.df = get.z.func(blood3.exp.control.mean,blood3.exp.control.sd,blood3.exp.asd)
blood3.sig.genes = gene.name.func(blood3.z.df)
blood3.sig.gene.sym = get.gene.symbol.func(blood3.sig.genes,blood3.info)
#####################################################################################################
# separate perturbed genes by sex
blood3.sex = gsub(".*, ","",blood3.pdat$source_name_ch1)
blood3.sex = blood3.sex[blood3.pdat$diagnosis == "Autism"]
blood3.male.genes = blood3.sig.gene.sym[blood3.sex=="Male"]
blood3.female.genes = blood3.sig.gene.sym[blood3.sex=="Female"]
rm(blood3.exp.control.mean,blood3.exp.control.sd,blood3.sig.genes)
# identify differentially expressed genes via traditional t-test
blood3.pvals = apply(blood3.exp,1,function(x){t.test(x~blood3.ed)$p.value})
blood3.sig = blood3.pvals[blood3.pvals<0.05]
blood3.sig.names = names(blood3.sig)
blood3.sig.names = unique(blood3.info[blood3.info$ID %in% blood3.sig.names,]$`Gene symbol`)
