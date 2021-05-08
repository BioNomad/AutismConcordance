library(GEOquery)
library(limma)
library(org.Hs.eg.db)
######################################################################################################
# first female data set
data_to_pull = read.table(file = "./science/autism_gene/blood/data_to_pull",
                          header = T,fill = T)
blood10 = getGEO(GEO= data_to_pull[,1][10],GSEMatrix = TRUE,AnnotGPL = T)
blood10.exp = as.data.frame(exprs(blood10$GSE123302_series_matrix.txt.gz))
blood10.pdat = as.data.frame(pData(blood10$GSE123302_series_matrix.txt.gz))
blood10.info = blood10$GSE123302_series_matrix.txt.gz@featureData@data
rm(blood10)
gene.names <- select(org.Hs.eg.db,keys=blood10.info$GB_ACC,
                     keytype="ACCNUM",columns=c("SYMBOL",
                                                 "ACCNUM"))
blood10.info$GENE_SYMBOL = gene.names$SYMBOL
blood10.exp = blood10.exp[,blood10.pdat$`diagnosis:ch1`=="ASD"|blood10.pdat$`diagnosis:ch1`=="TD"]
blood10.pdat = blood10.pdat[blood10.pdat$`diagnosis:ch1`=="ASD"|blood10.pdat$`diagnosis:ch1`=="TD",]
blood10.ed.tmp = blood10.pdat$`diagnosis:ch1`
blood10.ed.tmp = gsub("ASD",1,blood10.ed.tmp)
blood10.ed.tmp = gsub("TD",0,blood10.ed.tmp)
blood10.ed = factor(blood10.ed.tmp,levels = c(0,1))
#####################################################################################################
# identify genes perturbed per patient:
blood10.exp.asd = blood10.exp[blood10.ed ==1]
blood10.exp.control = blood10.exp[blood10.ed ==0]
blood10.exp.control.sd = apply(blood10.exp.control,1,sd)
blood10.exp.control.mean = apply(blood10.exp.control,1,mean)
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
    new.info = y[x[[i]],]$GENE_SYMBOL
    new.info = new.info[new.info>1]
    new.patient.list[[i]] = new.info
  }
  return(new.patient.list)
}
# apply the function to blood10 data 
blood10.z.df = get.z.func(blood10.exp.control.mean,blood10.exp.control.sd,blood10.exp.asd)
blood10.sig.genes = gene.name.func(blood10.z.df)
blood10.sig.gene.sym = get.gene.symbol.func(blood10.sig.genes,blood10.info)
######################################################################################################
rm(blood10.exp.control.mean,blood10.exp.control.sd,blood10.sig.genes)
