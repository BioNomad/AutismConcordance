library(DGCA)
source("./autism_concordance/scripts/common_perturbed.R")
##############################################################################################
# split by gender
# set up design matrices
# use spearman correlation to find which genes are highly correlated with 
# our unnamed genes
design_mat1.male = model.matrix(~ blood1.ed + 0)
colnames(design_mat1.male) = c("TD","ASD")
orf_id1.male = blood1.info[blood1.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res1.male = ddcorAll(inputMat = blood1.exp, design = design_mat1.male,
                      compare = c("TD", "ASD"),
                      adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                      corrType = "spearman",splitSet = orf_id1.male)
ddcor_res1.male = ddcor_res1.male[ddcor_res1.male$ASD_cor<(-.6)|ddcor_res1.male$ASD_cor>.6,]
ddcor_res1.male = ddcor_res1.male[ddcor_res1.male$pValDiff<0.05,]
ddcor_res1.male$Gene1=blood1.info[ddcor_res1.male$Gene1,]$`Gene symbol`
ddcor_res1.male$Gene2=blood1.info[ddcor_res1.male$Gene2,]$`Gene symbol`
ddcor_res1.male.high = ddcor_res1.male[ddcor_res1.male$ASD_cor<(-.7)|ddcor_res1.male$ASD_cor>.7,]


design_mat3.male = model.matrix(~ blood3.ed[gsub(".*, ","",blood3.pdat.select$source_name_ch1)=="Male"] + 0)
colnames(design_mat3.male) = c("TD","ASD")
orf_id3.male = blood3.info[blood3.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res3.male = ddcorAll(inputMat = blood3.exp[,gsub(".*, ","",blood3.pdat.select$source_name_ch1)=="Male"], design = design_mat3.male,
                      compare = c("TD", "ASD"),
                      adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                      corrType = "spearman",splitSet = orf_id3.male)
ddcor_res3.male = ddcor_res3.male[ddcor_res3.male$ASD_cor<(-0.6)|ddcor_res3.male$ASD_cor>(0.6),]
ddcor_res3.male = ddcor_res3.male[ddcor_res3.male$pValDiff<0.05,]
ddcor_res3.male$Gene1=blood3.info[ddcor_res3.male$Gene1,]$`Gene symbol`
ddcor_res3.male$Gene2=blood3.info[ddcor_res3.male$Gene2,]$`Gene symbol`
ddcor_res3.male.high = ddcor_res3.male[ddcor_res3.male$ASD_cor<(-.7)|ddcor_res3.male$ASD_cor>.7,]

design_mat8.male = model.matrix(~ blood8.ed[blood8.pdat$`gender:ch1`=="male"] + 0)
colnames(design_mat8.male) = c("TD","ASD")
orf_id8.male = blood8.info[blood8.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res8.male = ddcorAll(inputMat = blood8.exp[,blood8.pdat$`gender:ch1`=="male"], design = design_mat8.male,
                      compare = c("TD", "ASD"),
                      adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                      corrType = "spearman",splitSet = orf_id8.male)
ddcor_res8.male = ddcor_res8.male[ddcor_res8.male$ASD_cor<(-0.6)|ddcor_res8.male$ASD_cor>(0.6),]
ddcor_res8.male = ddcor_res8.male[ddcor_res8.male$pValDiff<0.05,]
ddcor_res8.male$Gene1=blood8.info[ddcor_res8.male$Gene1,]$`Gene symbol`
ddcor_res8.male$Gene2=blood8.info[ddcor_res8.male$Gene2,]$`Gene symbol`
ddcor_res8.male.high = ddcor_res8.male[ddcor_res8.male$ASD_cor<(-.7)|ddcor_res8.male$ASD_cor>.7,]

design_mat8.female = model.matrix(~ blood8.ed[blood8.pdat$`gender:ch1`=="female"] + 0)
colnames(design_mat8.female) = c("TD","ASD")
orf_id8.female = blood8.info[blood8.info$`Gene symbol`%in% female.orf.genes,]$ID
ddcor_res8.female = ddcorAll(inputMat = blood8.exp[,blood8.pdat$`gender:ch1`=="female"], design = design_mat8.female,
                           compare = c("TD", "ASD"),
                           adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                           corrType = "spearman",splitSet = orf_id8.female)
ddcor_res8.female = ddcor_res8.female[ddcor_res8.female$ASD_cor<(-0.6)|ddcor_res8.female$ASD_cor>(0.6),]
ddcor_res8.female = ddcor_res8.female[ddcor_res8.female$pValDiff<0.05,]
ddcor_res8.female$Gene1=blood8.info[ddcor_res8.female$Gene1,]$`Gene symbol`
ddcor_res8.female$Gene2=blood8.info[ddcor_res8.female$Gene2,]$`Gene symbol`
ddcor_res8.female.high = ddcor_res8.female[ddcor_res8.female$ASD_cor<(-.7)|ddcor_res8.female$ASD_cor>.7,]

design_mat9.male = model.matrix(~ blood9.ed + 0)
colnames(design_mat9.male) = c("TD","ASD")
orf_id9.male = blood9.info[blood9.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res9.male = ddcorAll(inputMat = blood9.exp, design = design_mat9.male,
                           compare = c("TD", "ASD"),
                           adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                           corrType = "spearman",splitSet = orf_id9.male)
ddcor_res9.male = ddcor_res9.male[ddcor_res9.male$ASD_cor<(-.6)|ddcor_res9.male$ASD_cor>.6,]
ddcor_res9.male = ddcor_res9.male[ddcor_res9.male$pValDiff<0.05,]
ddcor_res9.male$Gene1=blood9.info[ddcor_res9.male$Gene1,]$`Gene symbol`
ddcor_res9.male$Gene2=blood9.info[ddcor_res9.male$Gene2,]$`Gene symbol`
ddcor_res9.male.high = ddcor_res9.male[ddcor_res9.male$ASD_cor<(-.7)|ddcor_res9.male$ASD_cor>.7,]

design_mat9.1.male = model.matrix(~ blood9.1.ed[blood9.1.pdat$`gender:ch1`=="male"] + 0)
colnames(design_mat9.1.male) = c("TD","ASD")
orf_id9.1.male = blood9.1.info[blood9.1.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res9.1.male = ddcorAll(inputMat = blood9.1.exp[,blood9.1.pdat$`gender:ch1`=="male"], design = design_mat9.1.male,
                           compare = c("TD", "ASD"),
                           adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                           corrType = "spearman",splitSet = as.character(orf_id9.1.male))
ddcor_res9.1.male = ddcor_res9.1.male[ddcor_res9.1.male$ASD_cor<(-0.6)|ddcor_res9.1.male$ASD_cor>(0.6),]
ddcor_res9.1.male = ddcor_res9.1.male[ddcor_res9.1.male$pValDiff<0.05,]
ddcor_res9.1.male$Gene1=blood9.1.info[ddcor_res9.1.male$Gene1,]$`Gene symbol`
ddcor_res9.1.male$Gene2=blood9.1.info[ddcor_res9.1.male$Gene2,]$`Gene symbol`
ddcor_res9.1.male.high = ddcor_res9.1.male[ddcor_res9.1.male$ASD_cor<(-.7)|ddcor_res9.1.male$ASD_cor>.7,]

design_mat9.1.female = model.matrix(~ blood9.1.ed[blood9.1.pdat$`gender:ch1`=="female"] + 0)
colnames(design_mat9.1.female) = c("TD","ASD")
orf_id9.1.female = blood9.1.info[blood9.1.info$`Gene symbol`%in% female.orf.genes,]$ID
ddcor_res9.1.female = ddcorAll(inputMat = blood9.1.exp[,blood9.1.pdat$`gender:ch1`=="female"], design = design_mat9.1.female,
                             compare = c("TD", "ASD"),
                             adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                             corrType = "spearman",splitSet = as.character(orf_id9.1.female))
ddcor_res9.1.female = ddcor_res9.1.female[ddcor_res9.1.female$ASD_cor<(-0.6)|ddcor_res9.1.female$ASD_cor>(0.6),]
ddcor_res9.1.female = ddcor_res9.1.female[ddcor_res9.1.female$pValDiff<0.05,]
ddcor_res9.1.female$Gene1=blood9.1.info[ddcor_res9.1.female$Gene1,]$`Gene symbol`
ddcor_res9.1.female$Gene2=blood9.1.info[ddcor_res9.1.female$Gene2,]$`Gene symbol`
ddcor_res9.1.female.high = ddcor_res9.1.female[ddcor_res9.1.female$ASD_cor<(-.7)|ddcor_res9.1.female$ASD_cor>.7,]

design_mat11.male = model.matrix(~ blood11.ed + 0)
colnames(design_mat11.male) = c("TD","ASD")
orf_id11.male = blood11.info[blood11.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res11.male = ddcorAll(inputMat = blood11.exp, design = design_mat11.male,
                       compare = c("TD", "ASD"),
                       adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                       corrType = "spearman",splitSet = orf_id11.male)
ddcor_res11.male = ddcor_res11.male[ddcor_res11.male$ASD_cor<(-0.6)|ddcor_res11.male$ASD_cor>(0.6),]
ddcor_res11.male = ddcor_res11.male[ddcor_res11.male$pValDiff<0.05,]
ddcor_res11.male$Gene1=blood11.info[ddcor_res11.male$Gene1,]$`Gene symbol`
ddcor_res11.male$Gene2=blood11.info[ddcor_res11.male$Gene2,]$`Gene symbol`
ddcor_res11.male.high = ddcor_res11.male[ddcor_res11.male$ASD_cor<(-.7)|ddcor_res11.male$ASD_cor>.7,]

design_mat12.male = model.matrix(~ blood12.ed + 0)
colnames(design_mat12.male) = c("TD","ASD")
orf_id12.male = blood12.info[blood12.info$`Gene symbol`%in% male.orf.genes,]$ID
ddcor_res12.male = ddcorAll(inputMat = blood12.exp, design = design_mat12.male,
                       compare = c("TD", "ASD"),
                       adjust = "BH", heatmapPlot = FALSE, nPerm = 0, nPairs = 100,
                       corrType = "spearman",splitSet = orf_id12.male)
ddcor_res12.male = ddcor_res12.male[ddcor_res12.male$ASD_cor<(-0.6)|ddcor_res12.male$ASD_cor>(0.6),]
ddcor_res12.male = ddcor_res12.male[ddcor_res12.male$pValDiff<0.05,]
ddcor_res12.male$Gene1=blood12.info[ddcor_res12.male$Gene1,]$`Gene symbol`
ddcor_res12.male$Gene2=blood12.info[ddcor_res12.male$Gene2,]$`Gene symbol`
ddcor_res12.male.high = ddcor_res12.male[ddcor_res12.male$ASD_cor<(-.7)|ddcor_res12.male$ASD_cor>.7,]
########################################################################################################
# generate lists of differentially correlated genes per unnamed gene and gender
dgca.func = function(results,orf.genes){
  dgca.list = list()
  for (i in 1:length(orf.genes)){
    dgca.genes = results[results$Gene2==orf.genes[i],]$Gene1
    dgca.list[[i]] = dgca.genes
  }
  return(dgca.list)
}

blood1.dgca.male.high = dgca.func(ddcor_res1.male.high,male.orf.genes)
blood3.dgca.male.high = dgca.func(ddcor_res3.male.high,male.orf.genes)
blood8.dgca.male.high = dgca.func(ddcor_res8.male.high,male.orf.genes)
blood9.dgca.male.high = dgca.func(ddcor_res9.male.high,male.orf.genes)
blood9.1.dgca.male.high = dgca.func(ddcor_res9.1.male.high,male.orf.genes)
blood11.dgca.male.high = dgca.func(ddcor_res11.male.high,male.orf.genes)
blood12.dgca.male.high = dgca.func(ddcor_res12.male.high,male.orf.genes)

blood8.dgca.female.high = dgca.func(ddcor_res8.female.high,female.orf.genes)
blood9.1.dgca.female.high = dgca.func(ddcor_res9.1.female.high,female.orf.genes)
#########################################################################################
# corroborate these correlations by checking them against low grade glioma data
# and only return correlated genes that can be confirmed 
find.common.dgca.male.high = function(orfs){
  lst = list()
  for (i in 1:length(orfs)){
    tmp=table(c(blood1.dgca.male.high[[i]],
                blood3.dgca.male.high[[i]],
                blood8.dgca.male.high[[i]],
                blood9.dgca.male.high[[i]],
                blood9.1.dgca.male.high[[i]],
                blood11.dgca.male.high[[i]],
                blood12.dgca.male.high[[i]]
    ))
    tmp=names(tmp)
    file = paste("TGCAPanCancerLowGradeGlioma",orfs[[i]],".tsv",sep = "")
    coexp = read.csv(paste("./autism_concordance/data/",file,sep=""),
                     sep = "\t",header=T)
    coexp = coexp[coexp$p.Value<0.05,]
    coexp = coexp[coexp$Spearman.s.Correlation<(-.3)|coexp$Spearman.s.Correlation>(0.3),]
    tmp.select = tmp[tmp %in% coexp$Correlated.Gene]
    lst[[i]] = tmp.select
  }
  return(lst)
}
common.dgca.male.high.list1 = find.common.dgca.male.high(male.orf.genes)


find.common.dgca.female.high = function(orfs){
  lst = list()
  for (i in 1:length(orfs)){
    tmp=table(c(blood8.dgca.female.high[[i]],
                blood9.1.dgca.female.high[[i]]
    ))
    tmp = names(tmp)
    file = paste("TGCAPanCancerLowGradeGlioma",orfs[[i]],".tsv",sep = "")
    coexp = read.csv(paste("./autism_concordance/data/",file,sep=""),
                     sep = "\t",header=T)
    coexp = coexp[coexp$p.Value<0.05,]
    coexp = coexp[coexp$Spearman.s.Correlation<(-.3)|coexp$Spearman.s.Correlation>(0.3),]
    tmp.select = tmp[tmp %in% coexp$Correlated.Gene]
    lst[[i]] = tmp.select
  }
  return(lst)
}
common.dgca.female.high.list1 = find.common.dgca.female.high(female.orf.genes)
