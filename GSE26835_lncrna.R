gse26835 <- read.csv(file = '/Users/Weian/Documents/Lab/lncRNA/Data/GSE26835/gse26835_matrix_rma_done.csv',header = TRUE,row.names = 1)
lncrna <- read.csv(file = '/Users/Weian/Documents/Lab/lncRNA/lncRNA_radiation/hgu133av2_lncrna_list.csv')


# quantile
library(preprocessCore)
gse26835_qn <- normalize.quantiles(as.matrix(gse26835)) 
row.names(gse26835_qn) <- row.names(gse26835)
colnames(gse26835_qn) <- colnames(gse26835)

lnc_probe <- which(row.names(gse26835_qn) %in% lncrna$Probe)
gse26835_qn.lnc <- gse26835_qn[lnc_probe,]
gse26835_qn.nolnc <- gse26835_qn[-lnc_probe,]


Control <- gse26835_qn_lnc[,grep('0hr',colnames(gse26835_qn_lnc))]
Rad_2hr <- gse26835_qn_lnc[,grep('2hr',colnames(gse26835_qn_lnc))]
Rad_6hr <- gse26835_qn_lnc[,grep('6hr',colnames(gse26835_qn_lnc))]
all <- data.frame(Control,Rad_2hr,Rad_6hr)


#anova
group <- factor(rep(c('Non_IR','Rad_2hr','Rad_6hr'),each = 207))
for(i in 1:dim(all)[1]){
  anova <- aov(as.numeric(all[i,1:621])~group)
  all$anova_p[i] <- summary(anova)[[1]]["Pr(>F)"][[1]][1]  
}
 
#TukeyHSD
gse26835_HSD <- all[which(all$anova_p <0.00001), ]
for(i in 1:dim(gse12626_HSD)[1]){
  anova <- aov(as.numeric(gse26835_HSD[i,1:621])~group)
  T_HSD <- TukeyHSD(anova)
  gse26835_HSD$Rad_Non_2hr[i] <- T_HSD[[1]][1,4]
  gse26835_HSD$Rad_Non_6hr[i] <- T_HSD[[1]][2,4]
  gse26835_HSD$Rad_6hr_2hr[i] <- T_HSD[[1]][3,4]
}



# #pairwise t-test (wrong process)
# pairtest <- list()
# for(i in 1:dim(all)[1]){
#   pairtest[[i]] <- pairwise.t.test(as.numeric(all[i,1:621]),group, p.adjust.method = "none")$p.value
# }
# 
# pairtest_unlist <- unlist(pairtest)
# all$Non_2hr = pairtest_unlist[seq(1,length(pairtest_unlist),4)]
# all$Non_6hr = pairtest_unlist[seq(2,length(pairtest_unlist),4)]
# all$IR2hr_6hr = pairtest_unlist[seq(4,length(pairtest_unlist),4)]
# 
# all_anova <- all[which(all$anova_p < 0.01),] # first anova, then goes below...
# all_anova_non_2hr <- all_anova[which(all_anova$Non_2hr < 0.01),] # Non vs 2hr
# all_anova_non_6hr <- all_anova[which(all_anova$Non_6hr < 0.01),] # Non vs 6hr
# all_anova_2hr_6hr <- all_anova[which(all_anova$IR2hr_6hr < 0.01),] # 6hr vs 2hr

gse26835_probe_non_2hr <- row.names(gse26835_HSD[which(gse26835_HSD$Rad_Non_2hr < 0.001),])
gse26835_probe_non_6hr <- row.names(gse26835_HSD[which(gse26835_HSD$Rad_Non_6hr < 0.001),])
gse26835_probe_2hr_6hr <- row.names(gse26835_HSD[which(gse26835_HSD$Rad_6hr_2hr < 0.001),])


library(gplots)
probe_venn <- list(Non_2hr = gse26835_probe_non_2hr, Non_6hr = gse26835_probe_non_6hr, IR_2hr_6hr = gse26835_probe_2hr_6hr)
a <- venn(probe_venn, small = 1.0) 


#------------t-test(old version)----------
radio_0_2 <- data.frame(cbind(Control,Rad_2hr))
radio_0_6 <- data.frame(cbind(Control,Rad_6hr))
radio_0_2_6 <- data.frame(cbind(Control,Rad_2hr,Rad_6hr))
# control vs 2hr
for(i in 1:dim(radio_0_2)[1]){
  radio_0_2$P_value[i] <- t.test(radio_0_2[i,1:155],radio_0_2[i,156:310],paired = TRUE)$p.value
  radio_0_2$FoldChange_2_0_log2[i] <- log2(mean(as.numeric(radio_0_2[i,1:155])) / mean(as.numeric(radio_0_2[i,156:310])))
}

radio_0_2.p <- radio_0_2[which(radio_0_2$P_value<0.05),]
cat('number of probe sets(0,2hr): ',dim(radio_0_2.p)[1])
plot(x = radio_0_2$FoldChange_2_0_log2, y = -log10(radio_0_2$P_value), xlab ="log2 FoldChange", ylab = "-log10 P-value")
points(x = radio_0_2.p$FoldChange_2_0_log2, y = -log10(radio_0_2.p$P_value),col='red',pch=12)

# control vs 6hr
for(i in 1:dim(radio_0_6)[1]){
  radio_0_6$P_value[i] <- t.test(radio_0_6[i,1:155],radio_0_6[i,156:310],paired = TRUE)$p.value
  radio_0_6$FoldChange_6_0_log2[i] <- log2(mean(as.numeric(radio_0_6[i,1:155])) / mean(as.numeric(radio_0_6[i,156:310])))
}
radio_0_6$Gene <- annot$Gene.Symbol

radio_0_6.p <- radio_0_6[which(radio_0_6$P_value < 0.0001),]
cat('number of probe sets(0,6hr): ',dim(radio_0_6.p)[1])

library(gplots)
heatmap.2(as.matrix(radio_0_6.p[,1:310]),col=redgreen(75), scale="row", Colv = FALSE,
          key=TRUE, symkey=FALSE, density.info="none", 
          trace="none", cexRow=0.5)


library(rgl)
t_radio_0_6.p <- data.frame(t(radio_0_6.p[,1:310]))
t_radio_0_6.p$Type <-c(rep('red',155),rep('blue',155))
fit <- princomp(t_radio_0_6.p[,1:42])
plot3d(fit$scores[,1:3],col = t_radio_0_6.p$Type)

library(ggplot2)
pca <- prcomp(radio_0_2_6)
data_pca <- data.frame(pca$x[,1:3],type = c(rep('Non_IR',155),rep('IR_2hr',155),rep('IR_6hr',155)))


