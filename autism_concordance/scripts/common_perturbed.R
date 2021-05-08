source("./autism_concordance/scripts/blood1.R")
source("./autism_concordance/scripts/blood3.R")
source("./autism_concordance/scripts/blood8.R")
source("./autism_concordance/scripts/blood9.R")
source("./autism_concordance/scripts/blood91.R")
source("./autism_concordance/scripts/blood10.R")
source("./autism_concordance/scripts/blood11.R")
source("./autism_concordance/scripts/blood12.R")
library(RColorBrewer)
library(ggplot2)
library(sur)
##############################################################################
# find commonly perturbed genes between each study's perturbed gene pool
# also find commonly perturbed genes between males and females
common.perturbed = Reduce(intersect,list(
  unique(unlist(blood1.sig.gene.sym)),
  unique(unlist(blood3.sig.gene.sym)),
  unique(unlist(blood8.sig.gene.sym)),
  unique(unlist(blood9.sig.gene.sym)),
  unique(unlist(blood9.1.sig.gene.sym)),
  unique(unlist(blood11.sig.gene.sym)),
  unique(unlist(blood12.sig.gene.sym))
))
common.perturbed.male = Reduce(intersect,list(
  unique(unlist(blood1.sig.gene.sym)),
  unique(unlist(blood3.male.genes)),
  unique(unlist(blood8.male.genes)),
  unique(unlist(blood9.sig.gene.sym)),
  unique(unlist(blood9.1.male.genes)),
  unique(unlist(blood11.sig.gene.sym)),
  unique(unlist(blood12.sig.gene.sym))
))
common.perturbed.female = Reduce(intersect,list(
  unique(unlist(blood3.female.genes)),
  unique(unlist(blood8.female.genes)),
  unique(unlist(blood9.1.female.genes))
))
# find unnamed genes in each perturbed pool
orf.genes = common.perturbed[grepl("orf",common.perturbed)]
female.orf.genes = common.perturbed.female[grepl("orf",common.perturbed.female)]
male.orf.genes = common.perturbed.male[grepl("orf",common.perturbed.male)]
all.orf.genes = c(male.orf.genes,female.orf.genes)
#################################################################################
# determine which patients have these orf genes perturbed
get.common.func = function(common,lst){
  new.lst = list()
  for (i in 1:length(lst)){
    tmp = common[common %in% lst[[i]]]
    new.lst[[i]] = tmp
  }
  return(new.lst)
}
blood1.orf = get.common.func(all.orf.genes,blood1.sig.gene.sym)
blood3.orf = get.common.func(all.orf.genes,blood3.sig.gene.sym)
blood8.orf = get.common.func(all.orf.genes,blood8.sig.gene.sym)
blood9.orf = get.common.func(all.orf.genes,blood9.sig.gene.sym)
blood9.1.orf = get.common.func(all.orf.genes,blood9.1.sig.gene.sym)
blood11.orf = get.common.func(all.orf.genes,blood11.sig.gene.sym)
blood12.orf = get.common.func(all.orf.genes,blood12.sig.gene.sym)
#################################################################################
# identify the percentage of patients per stuyd with these unnammed genes
colors = brewer.pal(9,"Pastel1")
filled.out.perc.func = function(x){
  count = logical()
  for (i in 1:length(x)){
    filled.tmp = x[[i]]
    filled = filled.tmp > 1
    count = c(count,filled)
  }
  count.perc = percent.table(count)
  return(count.perc)
}

filled.out.perc.orf = c(
  percent.table(unlist(lapply(blood1.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood3.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood8.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood9.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood9.1.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood11.orf,function(x){!is.na(x[1])}))),
  percent.table(unlist(lapply(blood12.orf,function(x){!is.na(x[1])})))
)
filled.out.perc.orf = filled.out.perc.orf[names(filled.out.perc.orf) == "TRUE"]
orf.df = as.data.frame(cbind(
  c("GSE26415","GSE6575","GSE37772","GSE18123-GPL570","GSE18123-GPL6244","GSE42133","GSE25507"),
  round(filled.out.perc.orf,digits = 1)
))
colnames(orf.df) = c("StudyName", "Percent")

ggplot(orf.df, aes(y=Percent, x=StudyName)) + 
  geom_bar(stat="identity",fill=colors[1],width = .8)+
  xlab("")+
  geom_text(aes(label=paste(round(filled.out.perc.orf,digits = 1),"%",sep = "")), vjust=1.6, size=3.5)+
  ggtitle(paste("Percent of Patients with a Personally Perturbed",
                "Orf Gene per Study's Gene Pool",sep = "\n"))+
  theme(axis.text = element_text(size = 9,angle = 35,hjust = 1),plot.title = element_text(hjust = 0.5))











