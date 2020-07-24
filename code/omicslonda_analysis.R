library(devtools)
library(zoo)
library(lubridate)
library(ggplot2)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)
library(nlme)
library(imputeTS)
require(mgcv)
library(SummarizedExperiment)
library(grid)
library(devtools)
library(zoo)
library(lubridate)
library(ggplot2)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)
library(nlme)
library(imputeTS)
library(grid)
library(RColorBrewer)
library(dendextend)
library(circlize)
library(cluster)
library(fpc)
library(dendextend)
library(stats)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ggtern)
library(rmcorr)
library(pheatmap)
library(viridis)
library(nlme)
library(imputeTS)
require(mgcv)
library(SummarizedExperiment)
library(ggtern)
library(ggmap)
library("dplyr")
library("forcats")



rm(list=ls())
set.seed(72626)




###############################################
########### Load data and metadata     ########
###############################################
load("data/Multi_Omics_Seasonal.RData", envir = parent.frame(), verbose = FALSE)
ipop_metadata_all = merge(ls, sc, by.x = "SubjectID", by.y = "SubjectID")
rownames(ipop_metadata_all) = ipop_metadata_all$SampleID



#############################
######## OmicsLonDA  ########
#############################
source("Dropbox/OmicsLonDA_Bioconductor/OmicsLonDA_GAMM_05122020/R/omicslonda.R")
source("Dropbox/OmicsLonDA_Bioconductor/OmicsLonDA_GAMM_05122020/R/omicslondaHelper.R")
source("Dropbox/OmicsLonDA_Bioconductor/OmicsLonDA_GAMM_05122020/R/omicslondaMCPermutation.R")
source("Dropbox/OmicsLonDA_Bioconductor/OmicsLonDA_GAMM_05122020/R/omicslondaVisualization.R")


#######################################
########### Gut Microbiome     ########
#######################################
load("data/ipop_st_clr_abundant.RData")
ipop_gut_microbiome = ipop_st_clr_abundant
ipop_gut_microbiome = ipop_gut_microbiome[which(rownames(ipop_gut_microbiome) %in% sample_names), ]
ipop_gut_microbiome_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_gut_microbiome)), ]
colnames(ipop_gut_microbiome_metadata)[which(colnames(ipop_gut_microbiome_metadata) == "SubjectID")] = "Subject"
colnames(ipop_gut_microbiome_metadata)[which(colnames(ipop_gut_microbiome_metadata) == "IRIS")] = "Group"
colnames(ipop_gut_microbiome_metadata)[which(colnames(ipop_gut_microbiome_metadata) == "CollectionDate")] = "Date"
colnames(ipop_gut_microbiome_metadata)[which(colnames(ipop_gut_microbiome_metadata) == "CL4")] = "Status"
ipop_gut_microbiome_metadata$Date = as.Date(ipop_gut_microbiome_metadata$Date , format = '%m/%d/%y')
ipop_gut_microbiome_metadata$Time  = yday(ipop_gut_microbiome_metadata$Date)


dim(ipop_gut_microbiome_metadata)
dim(ipop_gut_microbiome)

## Reorder dataframes
ipop_gut_microbiome_metadata = ipop_gut_microbiome_metadata[rownames(ipop_gut_microbiome),]
dim(ipop_gut_microbiome_metadata)




se_ome_matrix = as.matrix(t(ipop_gut_microbiome))
se_metadata = DataFrame(ipop_gut_microbiome_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)

for(i in 1:ncol(ipop_gut_microbiome)){
  pfx = "OmicsLonDA_Seasonal_gut_microbiome_ccs_all_correctStatus_p100_05122020"
  omicslonda_test_object = omicslonda_se_object[i,]
  
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  
  if(res == "noFit"){
    next
  }
  
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  
  visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                         xlabel = "Day of the Year",
                         ylabel = "CLR Normalized Count", 
                         col = c("#fc8d62", "#8da0cb"),
                         prefix = pfx)
  
  
  
  if(length(res$start) > 0 ){
    
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days", ylabel = "CLR Normalized Count", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "CLR Normalized Count",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    saveRDS(res, file = sprintf("%s/Feature_%s_results_%s.rds",
                                prefix = pfx, 
                                text = res$details$feature[1], 
                                fit.method = "ccs"))
  }
  
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}




#########################################
############# LMM for microbiome  #######
#########################################
microbes_selected_omicslonda =  unique(do.call(rbind, microbiome)$feature)
ipop_gut_microbiome_selected = ipop_gut_microbiome[,microbes_selected_omicslonda]
ipop_gut_microbiome_metadata$month = as.factor(month(ipop_gut_microbiome_metadata$Date))
ipop_gut_microbiome_metadata$Group = as.factor(ipop_gut_microbiome_metadata$Group)
df_merge = merge(ipop_gut_microbiome_selected, ipop_gut_microbiome_metadata, by.x = 'row.names', by.y = "SampleID")
df_merge_exercise  = df_merge[-which(is.na(df_merge$totalmetminweek)),]


lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_exercise = NA)
for(i in colnames(ipop_gut_microbiome_selected)){
  print(i)
  ml <- as.formula( paste( i, "~ Group + month + totalmetminweek") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge_exercise, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_exercise = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_exercise))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_exercise = as.numeric(lmm_df$pval_exercise)

lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_exercise = p.adjust(lmm_df$pval_exercise, method = "fdr")
write.csv(lmm_df, "lmm_df_microbiome_exercise.csv", row.names = FALSE)





#######################################
########### Clinical Data     #########
#######################################
ipop_clinical = clinic.df
rownames(ipop_clinical) = ipop_clinical$SampleID
ipop_clinical = ipop_clinical[,-which(colnames(ipop_clinical) %in% c("SampleID", "SubjectID", "CollectionDate", "CL1", "CL2", "CL3", "CL4"))]
ipop_clinical = ipop_clinical[which(rownames(ipop_clinical) %in% sample_names), ]
ipop_clinical_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_clinical)), ]
colnames(ipop_clinical_metadata)[which(colnames(ipop_clinical_metadata) == "SubjectID")] = "Subject"
colnames(ipop_clinical_metadata)[which(colnames(ipop_clinical_metadata) == "IRIS")] = "Group"
colnames(ipop_clinical_metadata)[which(colnames(ipop_clinical_metadata) == "CollectionDate")] = "Date"
colnames(ipop_clinical_metadata)[which(colnames(ipop_clinical_metadata) == "CL4")] = "Status"
ipop_clinical_metadata$Date = as.Date(ipop_clinical_metadata$Date , format = '%m/%d/%y')
ipop_clinical_metadata$Time  = yday(ipop_clinical_metadata$Date)


####### Visualize Spline 
feature = "A1C"
x = as.data.frame(ipop_clinical[,c("SampleID", feature)])
df = merge(x, ipop_metadata_all,  by = "SampleID")

pdf("A1C.pdf", h = 5, w = 7)
ggplot(df, aes(Time, A1C)) + 
  geom_point(stat='identity', size = 1, show.legend = TRUE, alpha = 0.6, color = "black") +
  theme_bw() + 
  geom_smooth() +
  labs(y = "A1C level", x = "Day of the Year") +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size=15, face="bold"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"))
dev.off()



dim(ipop_clinical_metadata)
dim(ipop_clinical)

## Reorder dataframes
### TODO: Check for duplicates in names
ipop_clinical_metadata = ipop_clinical_metadata[rownames(ipop_clinical),]


# remove features that have missing values > 5% of ssample size
rm_feature = vector()
for(i in 1:ncol(ipop_clinical)){
  print(colnames(ipop_clinical)[i])
  
  x = ipop_clinical[,i]
  len = sum(is.na(x) | x == "")
  
  print(len)
  if(len>0.05*nrow(ipop_clinical)){
    rm_feature = c(rm_feature, i)
  }
}


ipop_clinical = ipop_clinical[,-rm_feature]
ipop_clinical = apply(ipop_clinical, 2, function(x) gsub("^$|^ $", NA, x))
ipop_clinical = apply(ipop_clinical, 2, function(x) gsub(">", "", x))
ipop_clinical = apply(ipop_clinical, 2, function(x) gsub("<", "", x))

temp = apply(ipop_clinical, 2, as.numeric)
rownames(temp) = rownames(ipop_clinical)
ipop_clinical = temp

se_ome_matrix = as.matrix(t(ipop_clinical))
se_metadata = DataFrame(ipop_clinical_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)


for(i in 1:ncol(ipop_clinical)){
  pfx = "OmicsLonDA_Seasonal_clinical_ccs_all_correctStatus_p100"
  omicslonda_test_object = omicslonda_se_object[i,]
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  
  if(res == "noFit"){
    next
  }
  
  if(length(res$start) > 0 ){
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days", ylabel = "Level", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                           xlabel = "Day of the Year",
                           ylabel = "Level", 
                           col = c("#fc8d62", "#8da0cb"),
                           prefix = pfx)
    
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "Level",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    save(res, file = sprintf("%s/Feature_%s_results_%s.RData",
                             prefix = pfx, 
                             text = res$details$feature[1], 
                             fit.method = "ccs"))
  }
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}






#########################################
############# LMM for clinical  ########
#########################################
selected_omicslonda =  unique(do.call(rbind, clinical)$feature)
ipop_ome_selected = ipop_clinical[,selected_omicslonda]
ipop_clinical_metadata$month = as.factor(month(ipop_clinical_metadata$Date))
ipop_clinical_metadata$Group = as.factor(ipop_clinical_metadata$Group)
df_merge = merge(ipop_ome_selected, ipop_clinical_metadata, by.x = 'row.names', by.y = "SampleID")
df_merge_exercise  = df_merge[-which(is.na(df_merge$totalmetminweek)),]

lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_exercise = NA)
for(i in colnames(ipop_ome_selected)){
  print(i)
  df_merge_exercise[,i] =  as.numeric(df_merge_exercise[,i])
  ml <- as.formula( paste( i, "~ Group + month + totalmetminweek") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge_exercise, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_exercise = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_exercise))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_exercise = as.numeric(lmm_df$pval_exercise)

lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_exercise = p.adjust(lmm_df$pval_exercise, method = "fdr")
write.csv(lmm_df, "lmm_df_clinical_exercise.csv", row.names = FALSE)








#######################################
############# Cytokines    ############
#######################################
ipop_cytokines = ck.df
rownames(ipop_cytokines) = ipop_cytokines$SampleID
ipop_cytokines = ipop_cytokines[,-which(colnames(ipop_cytokines) %in% c("SampleID", "SampleName", "Plate", "SubjectID", "CollectionDate", "CL1", "CL2", "CL3", "CL4"))]
ipop_cytokines = ipop_cytokines[which(rownames(ipop_cytokines) %in% sample_names), ]
ipop_cytokines_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_cytokines)), ]
colnames(ipop_cytokines_metadata)[which(colnames(ipop_cytokines_metadata) == "SubjectID")] = "Subject"
colnames(ipop_cytokines_metadata)[which(colnames(ipop_cytokines_metadata) == "IRIS")] = "Group"
colnames(ipop_cytokines_metadata)[which(colnames(ipop_cytokines_metadata) == "CollectionDate")] = "Date"
colnames(ipop_cytokines_metadata)[which(colnames(ipop_cytokines_metadata) == "CL4")] = "Status"
ipop_cytokines_metadata$Date = as.Date(ipop_cytokines_metadata$Date , format = '%m/%d/%y')
ipop_cytokines_metadata$Time  = yday(ipop_cytokines_metadata$Date)


dim(ipop_cytokines_metadata)
dim(ipop_cytokines)

## Reorder dataframes
ipop_cytokines_metadata = ipop_cytokines_metadata[rownames(ipop_cytokines),]

se_ome_matrix = as.matrix(t(ipop_cytokines))
se_metadata = DataFrame(ipop_cytokines_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)

for(i in 1:ncol(ipop_cytokines)){
  pfx = "OmicsLonDA_Seasonal_cytokines_ccs_all_correctStatus_p100"
  omicslonda_test_object = omicslonda_se_object[i,]
  
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  
  if(res == "noFit"){
    next
  }
  
  
  if(length(res$start) > 0 ){
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days", ylabel = "Normalized Intensity", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                           xlabel = "Day of the Year",
                           ylabel = "Normalized Intensity", 
                           col = c("#fc8d62", "#8da0cb"),
                           prefix = pfx)
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "Normalized Intensity",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    save(res, file = sprintf("%s/Feature_%s_results_%s.RData",
                             prefix = pfx, 
                             text = res$details$feature[1], 
                             fit.method = "ccs"))
  }
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}


####### Visualize Spline 
feature = "IP10"
x = as.data.frame(ipop_cytokines[,c("SampleID", feature)])
df = merge(x, ipop_metadata_all,  by = "SampleID")
dim(df)
df = df[which(df$IP10<300), ]
dim(df)

pdf("IP10.pdf", h = 5, w = 7)
ggplot(df, aes(Time, IP10)) + 
  geom_point(stat='identity', size = 1, show.legend = TRUE, alpha = 0.6, color = "black") +
  theme_bw() + 
  geom_smooth() +
  labs(y = "IP10 level", x = "Day of the Year") +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size=15, face="bold"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"))
dev.off()



#######################################
############# Proteomics     ##########
#######################################
ipop_proteomics = swathprot.PCR.df
rownames(ipop_proteomics) = ipop_proteomics$SampleID
ipop_proteomics = ipop_proteomics[,-which(colnames(ipop_proteomics) %in% c("SampleID", "SubjectID", "CollectionDate", "CL1", "CL2", "CL3", "CL4"))]
ipop_proteomics = ipop_proteomics[which(rownames(ipop_proteomics) %in% sample_names), ]
ipop_proteomics_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_proteomics)), ]
colnames(ipop_proteomics_metadata)[which(colnames(ipop_proteomics_metadata) == "SubjectID")] = "Subject"
colnames(ipop_proteomics_metadata)[which(colnames(ipop_proteomics_metadata) == "IRIS")] = "Group"
colnames(ipop_proteomics_metadata)[which(colnames(ipop_proteomics_metadata) == "CollectionDate")] = "Date"
colnames(ipop_proteomics_metadata)[which(colnames(ipop_proteomics_metadata) == "CL4")] = "Status"
ipop_proteomics_metadata$Date = as.Date(ipop_proteomics_metadata$Date , format = '%m/%d/%y')
ipop_proteomics_metadata$Time  = yday(ipop_proteomics_metadata$Date)


dim(ipop_proteomics_metadata)
dim(ipop_proteomics)

## Reorder dataframes
ipop_proteomics_metadata = ipop_proteomics_metadata[rownames(ipop_proteomics),]


se_ome_matrix = as.matrix(t(ipop_proteomics))
se_metadata = DataFrame(ipop_proteomics_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)
for(i in 246:ncol(ipop_proteomics)){
  pfx = "OmicsLonDA_Seasonal_protomics_ccs_all_correctStatus_p100"
  omicslonda_test_object = omicslonda_se_object[i,]
  
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  
  if(res == "noFit"){
    next
  }
  
  if(length(res$start) > 0 ){
    
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days", ylabel = "Normalized Intensity", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                           xlabel = "Day of the Year",
                           ylabel = "Normalized Intensity", 
                           col = c("#fc8d62", "#8da0cb"),
                           prefix = pfx)
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "Normalized Intensity",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    save(res, file = sprintf("%s/Feature_%s_results_%s.RData",
                             prefix = pfx, 
                             text = res$details$feature[1], 
                             fit.method = "ccs"))
  }

  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}







#########################################
############# LMM for Proteomics  ########
#########################################
df_merge = merge(ipop_proteomics, ipop_proteomics_metadata, by = 'row.names')
lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_status = NA)
for(i in colnames(ipop_proteomics)){
  print(i)
  ml <- as.formula( paste( i, "~ Group + Time + Status") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_status = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_status))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_status = as.numeric(lmm_df$pval_status)
hist(lmm_df$pval_group)
hist(lmm_df$pval_time)
hist(lmm_df$pval_status)
lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_status = p.adjust(lmm_df$pval_status, method = "fdr")

write.csv(lmm_df, "lmm_df_proteins.csv", row.names = FALSE)





#########################################
############# LMM for proteomics  ########
#########################################
selected_omicslonda =  unique(do.call(rbind, proteomics)$feature)
ipop_ome_selected = ipop_proteomics[,selected_omicslonda]
ipop_proteomics_metadata$month = as.factor(month(ipop_proteomics_metadata$Date))
ipop_proteomics_metadata$Group = as.factor(ipop_proteomics_metadata$Group)
df_merge = merge(ipop_ome_selected, ipop_proteomics_metadata, by.x = 'row.names', by.y = "SampleID")
df_merge_exercise  = df_merge[-which(is.na(df_merge$totalmetminweek)),]

lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_exercise = NA)
for(i in colnames(ipop_ome_selected)){
  print(i)
  df_merge_exercise[,i] =  as.numeric(df_merge_exercise[,i])
  
  ml <- as.formula( paste( i, "~ Group + month + totalmetminweek") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge_exercise, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_exercise = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_exercise))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_exercise = as.numeric(lmm_df$pval_exercise)

lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_exercise = p.adjust(lmm_df$pval_exercise, method = "fdr")
write.csv(lmm_df, "lmm_df_proteomics_exercise.csv", row.names = FALSE)




#######################################
########### Metabolomics     ##########
#######################################
ipop_metabolomics = metbcr.df
rownames(ipop_metabolomics) = ipop_metabolomics$SampleID
ipop_metabolomics = ipop_metabolomics[,-which(colnames(ipop_metabolomics) %in% c("SampleID", "SubjectID", "CollectionDate", "CL1", "CL2", "CL3", "CL4"))]
ipop_metabolomics = ipop_metabolomics[which(rownames(ipop_metabolomics) %in% sample_names), ]
ipop_metabolomics2 = log(ipop_metabolomics)

for(i in 1:ncol(ipop_metabolomics2)){
  x = metb.curated[which(metb.curated$Compounds_ID == colnames(ipop_metabolomics2)[i]), ]$HMDB
  x = gsub("\\|","_",x)
  colnames(ipop_metabolomics2)[i] = x
}

ipop_metabolomics2 = ipop_metabolomics2[,-which(colnames(ipop_metabolomics2)=="")]



ipop_metabolomics_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_metabolomics2)), ]
colnames(ipop_metabolomics_metadata)[which(colnames(ipop_metabolomics_metadata) == "SubjectID")] = "Subject"
colnames(ipop_metabolomics_metadata)[which(colnames(ipop_metabolomics_metadata) == "IRIS")] = "Group"
colnames(ipop_metabolomics_metadata)[which(colnames(ipop_metabolomics_metadata) == "CollectionDate")] = "Date"
colnames(ipop_metabolomics_metadata)[which(colnames(ipop_metabolomics_metadata) == "CL4")] = "Status"
ipop_metabolomics_metadata$Date = as.Date(ipop_metabolomics_metadata$Date , format = '%m/%d/%y')
ipop_metabolomics_metadata$Time  = yday(ipop_metabolomics_metadata$Date)


dim(ipop_metabolomics_metadata)
dim(ipop_metabolomics2)

## Rerder dataframes
ipop_metabolomics_metadata = ipop_metabolomics_metadata[rownames(ipop_metabolomics2),]

se_ome_matrix = as.matrix(t(ipop_metabolomics2))
se_metadata = DataFrame(ipop_metabolomics_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)

for(i in 1:ncol(ipop_metabolomics2)){
  pfx = "OmicsLonDA_Seasonal_metabolomics_ccs_all_correctStatus_p100"
  omicslonda_test_object = omicslonda_se_object[i,]
  
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  
  if(res == "noFit"){
    next
  }
  
  if(length(res$start) > 0 ){
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days of the Year", ylabel = "Normalized Intensity", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                           xlabel = "Day of the Year",
                           ylabel = "Normalized Intensity", 
                           col = c("#fc8d62", "#8da0cb"),
                           prefix = pfx)
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "Normalized Intensity",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    save(res, file = sprintf("%s/Feature_%s_results_%s.RData",
                             prefix = pfx, 
                             text = res$details$feature[1], 
                             fit.method = "ccs"))
  }
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}




#########################################
########## LMM for metabolomics  ########
#########################################
df_merge = merge(ipop_metabolomics, ipop_metabolomics_metadata, by = 'row.names')
lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_status = NA)
for(i in colnames(ipop_metabolomics)[-1]){
  print(i)
  ml <- as.formula( paste( i, "~ Group + Time + Status") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_status = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_status))
}


lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_status = as.numeric(lmm_df$pval_status)
hist(lmm_df$pval_group)
hist(lmm_df$pval_time)
hist(lmm_df$pval_status)
lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_status = p.adjust(lmm_df$pval_status, method = "fdr")


write.csv(lmm_df, "lmm_df_metabolites.csv", row.names = FALSE)




#########################################
################ LMM for RNAseq  ########
#########################################
selected_omicslonda =  unique(do.call(rbind, metabolomics)$feature)
ipop_ome_selected = ipop_metabolomics2[,selected_omicslonda]
ipop_metabolomics_metadata$month = as.factor(month(ipop_metabolomics_metadata$Date))
ipop_metabolomics_metadata$Group = as.factor(ipop_metabolomics_metadata$Group)
df_merge = merge(ipop_ome_selected, ipop_metabolomics_metadata, by.x = 'row.names', by.y = "SampleID")
df_merge_exercise  = df_merge[-which(is.na(df_merge$totalmetminweek)),]

lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_exercise = NA)
colnames(ipop_ome_selected) = gsub("-", "_", colnames(ipop_ome_selected))
colnames(df_merge_exercise) = gsub("-", "_", colnames(df_merge_exercise))

for(i in colnames(ipop_ome_selected)){
  print(i)
  df_merge_exercise[,i] =  as.numeric(df_merge_exercise[,i])
  ml <- as.formula( paste(i, "~ Group + month + totalmetminweek") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge_exercise, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_exercise = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_exercise))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_exercise = as.numeric(lmm_df$pval_exercise)

lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_exercise = p.adjust(lmm_df$pval_exercise, method = "fdr")
write.csv(lmm_df, "lmm_df_metabolites_exercise.csv", row.names = FALSE)








#######################################
################# RNAseq     ##########
#######################################
ipop_rnaseq = rnaseq.log.df
rownames(ipop_rnaseq) = ipop_rnaseq$SampleID
ipop_rnaseq = ipop_rnaseq[,-which(colnames(ipop_rnaseq) %in% c("SampleID", "SubjectID", "CollectionDate", "CL1", "CL2", "CL3", "CL4"))]
ipop_rnaseq = ipop_rnaseq[which(rownames(ipop_rnaseq) %in% sample_names), ]
ipop_rnaseq_metadata = ipop_metadata[which(ipop_metadata$SampleID %in% rownames(ipop_rnaseq)), ]
colnames(ipop_rnaseq_metadata)[which(colnames(ipop_rnaseq_metadata) == "SubjectID")] = "Subject"
colnames(ipop_rnaseq_metadata)[which(colnames(ipop_rnaseq_metadata) == "IRIS")] = "Group"
colnames(ipop_rnaseq_metadata)[which(colnames(ipop_rnaseq_metadata) == "CollectionDate")] = "Date"
colnames(ipop_rnaseq_metadata)[which(colnames(ipop_rnaseq_metadata) == "CL4")] = "Status"
ipop_rnaseq_metadata$Date = as.Date(ipop_rnaseq_metadata$Date , format = '%m/%d/%y')
ipop_rnaseq_metadata$Time  = yday(ipop_rnaseq_metadata$Date)


dim(ipop_rnaseq_metadata)
dim(ipop_rnaseq)

## Reorder dataframes
ipop_rnaseq_metadata = ipop_rnaseq_metadata[rownames(ipop_rnaseq),]
which(colnames(ipop_rnaseq) == "AR")
colnames(ipop_rnaseq)[1000]

## TODO: Select ccndidate genes before running OmicsLonDA
se_ome_matrix = as.matrix(t(ipop_rnaseq))
se_metadata = DataFrame(ipop_rnaseq_metadata)
omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
                                            colData = se_metadata)


for(i in 1:ncol(ipop_rnaseq)){
  pfx = "OmicsLonDA_Seasonal_rnaseq_ccs_all_correctStatus_p100"
  omicslonda_test_object = omicslonda_se_object[i,]
  
  points = seq(1, 365, length.out = 365)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 100,
                   fit.method = "ccs", points = points,
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH")
  if(res == "noFit"){
    next
  }
  
  if(length(res$start) > 0 ){
    if (!dir.exists(pfx)){
      dir.create(file.path(pfx))
    }
    
    visualizeFeature(se_object = omicslonda_test_object,
                     xlabel = "Days of the Year", ylabel = "Normalized Count", 
                     col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    visualizeFeatureSpline(se_object = omicslonda_test_object, omicslonda_object = res, fit.method = "ccs",
                           xlabel = "Day of the Year",
                           ylabel = "Normalized Count", 
                           col = c("#fc8d62", "#8da0cb"),
                           prefix = pfx)
    
    visualizeTestStatHistogram(omicslonda_object = res,
                               fit.method = "ccs", prefix = pfx)
    
    visualizeArea(omicslonda_object = res, fit.method = "ccs",
                  xlabel = "Day of the Year", 
                  ylabel = "Normalized Count",
                  col = c("#fc8d62", "#8da0cb"), prefix = pfx)
    
    save(res, file = sprintf("%s/Feature_%s_results_%s.RData",
                             prefix = pfx, 
                             text = res$details$feature[1], 
                             fit.method = "ccs"))
  }
  
  if (!dir.exists(pfx)){
    dir.create(file.path(pfx))
  }
  
  feature.summary = as.data.frame(do.call(cbind, res$details),
                                  stringsAsFactors = FALSE)
  
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv",
                                            prefix = pfx, text = res$details$feature[1],
                                            fit.method = "ccs"), row.names = FALSE)
}


#########################################
############# LMM for RNAseq  ###########
#########################################
selected_omicslonda =  unique(do.call(rbind, transcriptomics)$feature)
ipop_ome_selected = ipop_rnaseq[,selected_omicslonda]
ipop_rnaseq_metadata$month = as.factor(month(ipop_rnaseq_metadata$Date))
ipop_rnaseq_metadata$Group = as.factor(ipop_rnaseq_metadata$Group)
df_merge = merge(ipop_ome_selected, ipop_rnaseq_metadata, by.x = 'row.names', by.y = "SampleID")
df_merge_exercise  = df_merge[-which(is.na(df_merge$totalmetminweek)),]

lmm_df = data.frame(measure = NA, pval_group = NA, pval_time = NA, pval_exercise = NA)
colnames(ipop_ome_selected) = gsub("-", "_", colnames(ipop_ome_selected))
colnames(df_merge_exercise) = gsub("-", "_", colnames(df_merge_exercise))

for(i in colnames(ipop_ome_selected)){
  print(i)
  df_merge_exercise[,i] =  as.numeric(df_merge_exercise[,i])
  ml <- as.formula( paste(i, "~ Group + month + totalmetminweek") )
  m1.nlme = lme(ml, random = ~ 1|Subject, data = df_merge_exercise, na.action = na.exclude)
  m1.nlme_anova = anova(m1.nlme)
  pval_group = m1.nlme_anova$`p-value`[2]
  pval_time = m1.nlme_anova$`p-value`[3]
  pval_exercise = m1.nlme_anova$`p-value`[4]
  lmm_df = rbind(lmm_df, c(i, pval_group, pval_time, pval_exercise))
}

lmm_df = lmm_df[-1,]
lmm_df$pval_group = as.numeric(lmm_df$pval_group)
lmm_df$pval_time = as.numeric(lmm_df$pval_time)
lmm_df$pval_exercise = as.numeric(lmm_df$pval_exercise)

lmm_df$padj_group = p.adjust(lmm_df$pval_group, method = "fdr")
lmm_df$padj_time = p.adjust(lmm_df$pval_time, method = "fdr")
lmm_df$padj_exercise = p.adjust(lmm_df$pval_exercise, method = "fdr")
write.csv(lmm_df, "lmm_df_rnaseq_exercise.csv", row.names = FALSE)






####### Visualize Spline 
feature = "PER1"
x = as.data.frame(ipop_rnaseq[,c(feature)])
rownames(x) = rownames(ipop_rnaseq)
colnames(x) = feature

df = merge(x, ipop_metadata_all,  by.x = "row.names", by.y = "SampleID")
dim(df)
df = df[which(df$PER1 != 0),]
dim(df)


pdf("PER1.pdf", h = 5, w = 7)
ggplot(df, aes(Time, PER1)) + 
  geom_point(stat='identity', size = 1, show.legend = TRUE, alpha = 0.6, color = "black") +
  #facet_grid(rows = vars(variable), scales = "free_y", labeller = labeller(variable = variable.labs)) +
  theme_bw() + 
  geom_smooth() +
  labs(y = "PER1 Expression", x = "Day of the Year") +
  #ggtitle("Weather parameters") + 
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        legend.text = element_text(size=15, face="bold"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"))
dev.off()



###############################################
#### Combine csv files into RData or rds
###############################################
integrateFile = function(dir = NULL, omicType = NULL) {
  file_list =  list.files(path = dir, pattern="*.csv") 
  ome = list()
  for (i in file_list){
    x = paste(dir, i, sep = "")
    print(x)
    feature = strsplit(i, ".csv")[[1]][1]
    print(feature)
    df = read.csv(x)
    ome[[feature]] = df
  }
  return(ome)
}

clinical = integrateFile(dir = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/clinical/", omicType = "clinical")
transcriptomics = integrateFile(dir = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/genes/", omicType = "transcriptomics")
proteomics = integrateFile(dir = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/proteins/", omicType = "proteomics")
metabolomics = integrateFile(dir = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/metabolites/", omicType = "metabolomics")
microbiome = integrateFile(dir = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/microbes/", omicType = "microbiome")



OmicsLonDA_IR_IS = list(microbiome = microbiome, transcriptomics = transcriptomics, proteomics = proteomics, metabolomics = metabolomics, clinical = clinical)
saveRDS(OmicsLonDA_IR_IS, file = "OmicsLonDA_IR_IS.rds")



##############################################
#### Cluster based on actual testStat value
#############################################
extract_sig_omicslonda = function(omics_path = "NULL"){
  tmp = list.files(path = omics_path, pattern="Feature.+csv")
  myfiles = lapply(paste(omics_path,tmp, sep = "/"), read.csv)
  
  ome = data.frame()
  for(i in 1:length(myfiles)){
    print(i)
    feature = as.character(unique(myfiles[[i]]$feature))
    print(feature)
    padjust = myfiles[[i]]$adjusted.pvalue
    tm = padjust
    tm[padjust <= 0.05/2] = 1
    tm[padjust > 0.05/2] = 0
    tm = tm*myfiles[[i]]$testStat
    ome = rbind(ome, tm)
    rownames(ome)[i] = feature
  }
  
  colnames(ome) = myfiles[[1]]$points
  ome_filtered = ome[apply(ome, 1, function(x) !all(x==0)),]
  return(ome_filtered)
}

genes = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/genes")
proteins = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/proteins")
microbiome = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/microbes/")
metabolites = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/metabolites")
clinical = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/clinical/")
#cytokines = extract_sig_omicslonda(omics_path = "OmicsLonDA_analysis/omicslonda_sig_features_09032019/cytokines/") # 67 features

dim(genes)
dim(proteins)
dim(microbiome)
dim(metabolites)
dim(clinical)




#### Visualize as Heatmap
## TODO: Need to add zero column to microbiome to match zeros
microbiome$`354` = 0
all_omes = rbind(proteins, metabolites, microbiome, genes, clinical)
time = colnames(all_omes) #seq(1:ncol(z_all_vis_filtered)) 
annot_column = data.frame(Time = time)
group = c(rep("Proteins", nrow(proteins)), 
          rep("Metabolites", nrow(metabolites)), 
          rep("Gut_Microbes", nrow(microbiome)), 
          rep("Genes", nrow(genes)),
          rep("Clinical", nrow(clinical))) 
annot_row =  data.frame(Group = as.factor(group))
rownames(annot_row) = rownames(all_omes)

indx_wZero = unname(which(apply(all_omes, 1, function(x) any(x==0))))


## woZeros
all_omes_woZeros = all_omes[-indx_wZero, ]
group_woZeros = group[-indx_wZero]
annot_row_woZeros = data.frame(Group = as.factor(group_woZeros))
all_omes = all_omes_woZeros
annot_row = annot_row_woZeros 
rownames(annot_row) = rownames(all_omes)


# wZeros
all_omes_wZeros = all_omes[indx_wZero, ]
group_wZeros = group[indx_wZero]
annot_row_wZeros = data.frame(Group = as.factor(group_wZeros))
all_omes = all_omes_wZeros
annot_row = annot_row_wZeros 
rownames(annot_row) = rownames(all_omes)



paletteLength = 11
myColor = colorRampPalette(c("#8da0cb", "white", "#fc8d62"))(paletteLength)
myBreaks <- c(seq(min(all_omes), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(all_omes)/paletteLength, max(all_omes), length.out=floor(paletteLength/2)))



pdf("iPOP_OmicsLonDA_PartOfYear_05142020.pdf", width = 10, height = 10)
pheatmap(
  mat               = all_omes,
  color             = myColor,
  breaks = myBreaks, 
  legend = TRUE,
  cluster_cols = FALSE, 
  cluster_rows = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_row    = annot_row,
  annotation_colors = list(Group = c(Proteins = "#E76BF3", Metabolites = "lightskyblue", Gut_Microbes = "#00BF7D", 
                                     Genes = "#A3A500", Clinical  = "#F8766D"),                          
                           Time = c("white", "gray")),
  drop_levels       = TRUE,
  fontsize          = 6,
  angle_col         = 45,
  fontsize_col      = 1, 
  fontsize_row      = 8, 
  border_color      = "black",
  legend.cex = 19,
  main              = ""
)
grid.text("Day of the Year", y=0.01, gp=gpar(fontsize=10))
dev.off()


### Visualize number of significant molecules/microbes
sig_df_wZeros = as.data.frame(table(group_wZeros))
x = data.frame("clinical",0)
names(x) = c("group_wZeros", "Freq")
sig_df_wZeros = rbind(sig_df_wZeros, x)

## Without Zeros
sig_df_woZeros = as.data.frame(table(group_woZeros))
x = data.frame("gut_microbiome", 0)
names(x) = c("group_woZeros", "Freq")
sig_df_woZeros = rbind(sig_df_woZeros, x)


colnames(sig_df_wZeros) = c("omic", "Part of the year difference")
colnames(sig_df_woZeros) = c("omic", "All year difference")
w = merge(sig_df_woZeros, sig_df_wZeros, all = TRUE)
q = melt(w)
q$omic = as.character(q$omic)
q$omic[which(q$omic == "gut_microbiome")] = "Gut_Microbes"
q$omic[which(q$omic == "clinical")] = "Clinical"
q$omic[which(q$omic == "genes")] = "Genes"
q$omic[which(q$omic == "metabolites")] = "Metabolites"
q$omic[which(q$omic == "proteins")] = "Proteins"

pdf("ipop_seasonal_number_sig_features_iris_05122020.pdf", h = 4, w = 7)
ggplot(data = q, aes(x = omic, y = value, color = omic, fill = omic)) +
  geom_bar(colour="black", stat="identity") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.75, show.legend = FALSE) + 
  guides(fill=FALSE) + 
  theme_bw() +
  ylim(0,90)+
  labs(y = "# of Features", x = "") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=13, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=13, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, color = "black", face = "bold")) +
  facet_wrap(facets = vars(variable))
dev.off()
