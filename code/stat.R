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


####################################
###### Quantitative all iPOP   #####
#################################### 
#### Impute missing data for BMI
ipop_metadata_all$BMI = na_interpolation(ipop_metadata_all$BMI, option ="linear")

####### Remove samples without collection dates
remov = which(ipop_metadata_all$CollectionDate == "")
ipop_metadata_all = ipop_metadata_all[-remov, ]
ipop_metadata_all$CollectionDate = as.Date(ipop_metadata_all$CollectionDate , format = '%m/%d/%y')
ipop_metadata_all$Time  = yday(ipop_metadata_all$CollectionDate)
ipop_metadata_all$month  = month(ipop_metadata_all$CollectionDate)



### Add physical activities
exercise = read.csv("data/ipaq.csv")
exercise$SampleID = as.character(exercise$SampleID)
exercise_subset = exercise[, c("SampleID", "totalmetminweek")]
ipop_metadata_all = merge(ipop_metadata_all, exercise_subset, by = "SampleID", all.x = TRUE)


pdf("ipop_seasonal_timepoints_distribution_all_ipop.pdf", h = 10, w = 10)
ggplot(ipop_metadata_all, aes(x = SubjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "#fc8d62", "#8da0cb"),  breaks=c("IR", "IS")) +
  scale_x_discrete(name ="Subject")+ 
  geom_point(size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Day of the year", breaks = round(seq(0, max(ipop_metadata_all$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()


pdf(file = "iPOP_all_samples_histogram.pdf", h = 5, w = 7)
ggplot(ipop_metadata_all, aes(x=Time)) +
  geom_histogram(bins = 12, breaks=seq(0, 365, by=30), colour='red') +
  scale_colour_manual(values = c("#fc8d62", "#8da0cb"), breaks = c("IR", "IS")) + #scale_color_brewer(palette="Dark2") +
  #scale_fill_manual(values = c("#fc8d62", "#8da0cb"), breaks = c("IR", "IS")) +
  theme_minimal()+theme_classic()+theme(legend.position="top") +
  labs(y = "Frequency", x = "Day of the Year") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0,
                                   hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0,
                                   hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0,
                                    hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90,
                                    hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"),
        legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") + scale_x_continuous(breaks = waiver())
dev.off()


####################################
######### Quantitative IR/IS   #####
####################################
remov = which(ipop_metadata_all$IRIS == "Unknown" | is.na(ipop_metadata_all$IRIS))
ipop_metadata = ipop_metadata_all[-remov, ]
sample_names = ipop_metadata$SampleID

length(unique(ipop_metadata$SubjectID))

table(ipop_metadata[, c("IRIS")])
table(ipop_metadata[, c("Ethnicity")])
table(ipop_metadata[, c("Gender")])


length(unique(ipop_metadata$SubjectID))
table(ipop_metadata[, c("IRIS")])

x = ipop_metadata[, c("SubjectID", "IRIS", "Gender", "Adj.age", "BMI", "SSPG", "Ethnicity")]
x$SSPG = as.numeric(x$SSPG)
x = unique(x)
table(x[, c("IRIS")])
table(x[, c("Gender")])
table(x[, c("Ethnicity")])

ir = subset(x, IRIS == "IR")
is = subset(x, IRIS == "IS")

table(ir[, c("Gender")])
table(is[, c("Gender")])

table(ir[, c("Ethnicity")])
table(is[, c("Ethnicity")])

range(ir$Adj.age)
mean(ir$Adj.age)
sd(ir$Adj.age)

range(ir$BMI)
mean(ir$BMI)
sd(ir$BMI)

range(ir$SSPG)
mean(ir$SSPG)
sd(ir$SSPG)

range(is$SSPG)
range(is$Adj.age)
mean(is$Adj.age)
sd(is$Adj.age)

range(is$BMI, na.rm = TRUE)
mean(is$BMI, na.rm = TRUE)
sd(is$BMI, na.rm = TRUE)

mean(is$SSPG)
sd(is$SSPG)






##### Visualization IR/IS cohort
pdf(file = "iPOP_IRIS_samples_histogram.pdf", h = 7, w = 10)
ggplot(ipop_metadata, aes(x=Time, color=IRIS, fill = IRIS)) +
  geom_histogram(position="dodge", bins = 12, breaks=seq(0, 365, by=30)) +
  scale_colour_manual(values = c("#fc8d62", "#8da0cb"), breaks = c("IR", "IS")) +
  scale_fill_manual(values = c("#fc8d62", "#8da0cb"), breaks = c("IR", "IS")) +
  theme_minimal()+theme_classic()+theme(legend.position="top") +
  labs(y = "Frequency", x = "Day of the Year") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0,
                                   hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0,
                                   hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0,
                                    hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90,
                                    hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"),
        legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") + scale_x_continuous(breaks = waiver())
dev.off()



########### IR vs IS Statistics 
wilcox.test(ir$Adj.age, is$Adj.age)
wilcox.test(ir$BMI, is$BMI)
wilcox.test(ir$SSPG, is$SSPG)
tbl = table(x$Gender, x$IRIS)
chisq.test(tbl)
tbl = table(x$Ethnicity, x$IRIS)
chisq.test(tbl)