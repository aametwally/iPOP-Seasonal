library(Mfuzz)
library(ggplot2)
library(reshape2)


set.seed(837383)


###### ELbow method
data=read.table("Coeffs_mFuzz_Cyt_Met_RNA_Pro.txt",header=T,sep="\t")
dim(data)
data = data[,c(2:9)]
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=2)


data=read.table("GAM_CoeffsT_Seas_data_PRO_sig_mFuzzy.txt",header=F,sep="\t")
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
       + ylab="Within groups sum of squares",main="Assessing the Optimal Number of Clusters with the Elbow Method",pch=20, cex=2)



## Fuzzy c-means clustring
a = read.table("Coeffs_mFuzz_Cyt_Met_RNA_Pro.txt", header = TRUE, sep = "\t")
idx = which(!duplicated(a$Genes))
a_noDuplicates = a[idx,]
rownames(a_noDuplicates) = a_noDuplicates$Genes
a_noDuplicates = a_noDuplicates[,-1]
object = new("ExpressionSet", exprs = as.matrix(a_noDuplicates))
timeset.s = standardise(object)
cl = mfuzz(timeset.s, c=2, m=mestimate(timeset.s))


### Draw clusters with map to 0-365
clusters = as.data.frame(exprs(timeset.s)) #a_noDuplicates
clusters$cluster = cl$cluster
clusters$membership =  apply(cl$membership, 1, max) 
clusters$feature = rownames(clusters)
clust_membership = clusters[, c("feature", "membership")]


clust_1 = subset(clusters, cluster == 1)
x = clust_1[,-c(9,10)]
x_melt = melt(x)
x_melt$variable = as.numeric(gsub("t", "", x_melt$variable))
x_melt$variable = linMap(x_melt$variable, 0, 365)
clust = merge(x_melt, clust_membership, by.x = "feature", by.y = "feature")


jpeg("omics_pattern_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
ggplot(clust, aes(variable, value, colour = membership, group = feature)) + 
  theme_bw() + geom_point(size=1, alpha=0.4) + geom_line(linetype="solid", size=1, alpha=0.4) + 
  ggtitle("Pattern 1") + labs(y = "Coeff", x = "Day of the Year", fill = "Membership") +
  scale_colour_gradient(low = "white", high = "#1b9e77") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
        legend.text=element_text(size=8, face="plain"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(legend.position="right") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
dev.off()



clust_2 = subset(clusters, cluster == 2)
x = clust_2[,-c(9,10)]
x_melt = melt(x)
x_melt$variable = as.numeric(gsub("t", "", x_melt$variable))
x_melt$variable = linMap(x_melt$variable, 0, 365)
clust = merge(x_melt, clust_membership, by.x = "feature", by.y = "feature")


jpeg("omics_pattern_2.jpg", res = 300, height = 10, width = 15, units = 'cm')
ggplot(clust, aes(variable, value, colour = membership, group = feature)) + 
  theme_bw() + geom_point(size=1, alpha=0.4) + geom_line(linetype="solid", size=1, alpha=0.4) + 
  ggtitle("Pattern 2") + labs(y = "Coeff", x = "Day of the Year", fill = "Membership") +
  scale_colour_gradient(low = "white", high = "#d95f02") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
        legend.text=element_text(size=8, face="plain"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(legend.position="right") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
dev.off()



###### Correlation between seasons
a_noDuplicates$cluster = cl$cluster
a_noDuplicates_orderd = a_noDuplicates[order(a_noDuplicates$cluster),]
a_cor = cor(t(a_noDuplicates_orderd[, -ncol(a_noDuplicates_orderd)]))

melted_cormat <- melt(a_cor)
head(melted_cormat)

jpeg("ipop_patterns_correlation_03262019_new.jpg", res = 300, height = 75, width = 75, units = 'cm')
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x = element_text(colour="black",size=4,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=4,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  scale_fill_gradient2(low = "black", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")
dev.off()






### Create spline for each class
## Map function to convert the range from 1-8 to 0-365
linMap <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}


clust_1 = subset(a_noDuplicates_orderd, cluster == 1)
x = clust_1[,-9]
x_melt = melt(x)
x_melt$variable = as.numeric(gsub("t", "", x_melt$variable))
x_melt$variable = linMap(x_melt$variable, 0, 365)
loess_model = loess(value ~ variable, data = x_melt)
points = seq(min(x_melt$variable), max(x_melt$variable))
est.clust_1 = predict(loess_model, data.frame(variable = points), se = TRUE)
dd.est.clust_1 = data.frame(Time = points, Coeff = est.clust_1$fit, Group = "Pattern_1")




clust_2 = subset(a_noDuplicates_orderd, cluster == 2)
x = clust_2[,-9]
x_melt = melt(x)
x_melt$variable = as.numeric(gsub("t", "", x_melt$variable))
x_melt$variable = linMap(x_melt$variable, 0, 365)
loess_model = loess(value ~ variable, data = x_melt)
points = seq(min(x_melt$variable), max(x_melt$variable))
est.clust_2 = predict(loess_model, data.frame(variable = points), se = TRUE)
dd.est.clust_2 = data.frame(Time = points, Coeff = est.clust_2$fit, Group = "Pattern_2")


## Combine Datasets together
dd.est.clust = rbind(dd.est.clust_1, dd.est.clust_2)

jpeg("omics_seasons.jpg", res = 300, height = 10, width = 15, units = 'cm')
ggplot(dd.est.clust, aes(Time, Coeff, colour = Group)) + 
  theme_bw() + geom_point(size=1, alpha=0.9) + geom_line(linetype="solid", size=2, alpha=0.9) + 
  ggtitle("Omics Seasons") + labs(y = "Coeff", x = "Day of the Year") +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), 
                      breaks = c("Pattern_1", "Pattern_2", "Pattern_3")) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
dev.off()