library("mgcv")
library(Mfuzz)
library(ggplot2)
library(reshape2)
library(lubridate)
set.seed(747432)

#Input Data file containg features 
genes = read.table("Input_Data.txt",header = TRUE)
gene.names = names(genes)
genes = data.frame(genes)
genes = data.frame(sapply(genes, as.numeric))


#Input Annotation file describing Data file
ids = t(read.table("Input_colData.txt",header = FALSE))
colnames(ids) = ids[1,]
ids = data.frame(ids[-1,])
ids$Date = as.POSIXct(as.character(ids$Date), format="%e-%b-%Y")
ids$Time = yday(ids$Date)

ctrl <- lmeControl(opt='optim')
data = cbind(ids, genes)

p.vals = c()
coeffs = c()

# Build GAM models with mixed-effect
# The mixed-effect part removes the individual mean from the observed expression
# GAM models the average deviation from the mean of each subject
# it corresponds to the general trend of expression in the entire population
for (gene in gene.names) {
  png(paste0("plots/",gene,".png"))
  mod <- gamm(as.formula(paste(gene, " ~ IRIS + BMI + s(Time, bs = \"cc\")")),
              data = data, 
              method = "REML",
              control=ctrl,
              random = list(SubjectID = ~ 1),
              knots = list(TimeOfYear = c(0, 366)))
  
  plot(mod$gam, shade=T,shade.col="#9a9a00",xlab = "Day of the year", ylab = paste(gene),cex.lab=1.5,font.lab=2, cex.axis=1.5,font=2, main=paste(gene), font.main=2,cex.main=2.5)
  
  p.vals = rbind(p.vals, summary(mod$gam)$s.table[4])
  coeffsT = rbind(coeffsT, mod$gam$coefficients)
  
  
  dev.off()
}

p.vals_all = cbind(p.vals,gene.names)
coeffs_all = cbind(coeffs,gene.names)

write.table(p.vals_all, "Output_p.vals_GAMM.txt",sep="\t",quote=F)
write.table(coeffs_all, "Output_coeffs_GAMM.txt",sep="\t",quote=F)