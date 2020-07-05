library(limma)
datExpr<-read.csv('B26_CountMatrix_2509_Combined_sub289.csv',header=T,sep=',',row.names=1)
meta<-read.csv('meta289.csv',header=T,sep=',')
meta2<-meta[meta$missing==1,]
datExpr2<-datExpr[,meta$missing==1]

log2cpm<-voom(datExpr2)$E
nbp<-length(which(meta2$Diagnosis=='BP'))
nscz<-length(which(meta2$Diagnosis=='SCZ'))
nctl<-length(which(meta2$Diagnosis=='CTL'))

##gene filtering
thre<-apply(log2cpm,1,function(x){length(which(x[meta2$Diagnosis=='BP']>0.1))>nbp*0.75|length(which(x[meta2$Diagnosis=='SCZ']>0.1))>nscz*0.75|length(which(x[meta2$Diagnosis=='CTL']>0.1))>nctl*0.75})
count.fgene<-datExpr2[thre,]
log2cpm.fgene<-log2cpm[thre,]

##PCA plot
library(ggfortify)
library(gridExtra)
meta2$RIN<-as.numeric(meta2$RIN)
meta2$PMI<-as.numeric(meta2$PMI)
pc<-prcomp(t(log2cpm.fgene))
g1<-autoplot(pc, data = meta2, colour = 'Diagnosis')
g2<-autoplot(pc, data = meta2, colour = 'Batch')
g3<-autoplot(pc, data = meta2, colour = 'AgeDeath')
g4<-autoplot(pc, data = meta2, colour = 'Sex')
g5<-autoplot(pc, data = meta2, colour = 'Ethnicity')
g6<-autoplot(pc, data = meta2, colour = 'Collection')
g7<-autoplot(pc, data = meta2, colour = 'RIN')
g8<-autoplot(pc, data = meta2, colour = 'Hemisphere')
g9<-autoplot(pc, data = meta2, colour = 'PMI')
g10<-autoplot(pc, data = meta2, colour = 'BrainWeight')

pdf("log2cpm.fgene.pca.biological.pdf",width=10,height=8)
grid.arrange(g1,g3,g4,g5,g8,g10,nrow=2)
dev.off()


pdf("log2cpm.fgene.pca.technical.pdf",width=10,height=8)
grid.arrange(g1,g2,g6,g7,g5,g9,nrow=2)
dev.off()


#sample QC
library(WGCNA)
normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
#normadj <- (0.5+0.5*cor(log2cpm.fgene))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 3
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
log2cpm.fgene.fsample <- log2cpm.fgene[,!outliers]
count.fgene.fsample <- count.fgene[,!outliers]
mimeta<-meta2[!outliers,]
mimeta<-read.csv('mRNA/mimeta20200603.csv',header=T,sep=',')
save(log2cpm.fgene.fsample,count.fgene.fsample,mimeta,file='mirna289_fgene_fsample_sd3.RData')


##PCA after sample filtering
library(ggfortify)
library(gridExtra)
pc<-prcomp(t(log2cpm.fgene.fsample))
g1<-autoplot(pc, data = mimeta, colour = 'Diagnosis')
g2<-autoplot(pc, data = mimeta, colour = 'Batch')
g3<-autoplot(pc, data = mimeta, colour = 'AgeDeath')
g4<-autoplot(pc, data = mimeta, colour = 'Sex')
g5<-autoplot(pc, data = mimeta, colour = 'Ethnicity')
g6<-autoplot(pc, data = mimeta, colour = 'Collection')
g7<-autoplot(pc, data = mimeta, colour = 'RIN')
g8<-autoplot(pc, data = mimeta, colour = 'Hemisphere')
g9<-autoplot(pc, data = mimeta, colour = 'PMI')
g10<-autoplot(pc, data = mimeta, colour = 'BrainWeight')

pdf("log2cpm.fgene.fsample.pca.biological.pdf",width=10,height=8)
grid.arrange(g1,g3,g4,g5,g8,g10,nrow=2)
dev.off()


pdf("log2cpm.fgene.fsample.pca.technical.pdf",width=10,height=8)
grid.arrange(g1,g2,g6,g7,g5,g9,nrow=2)
dev.off()

##correlation between variables
mimeta$Diagnosis<-factor(mimeta$Diagnosis,levels=c('CTL','BP','SCZ'))
mimeta$Batch<-factor(mimeta$Batch)
mimeta$Hemisphere<-factor(mimeta$Hemisphere)
mimeta$Sex<-factor(mimeta$Sex)
mimeta$Ethnicity<-factor(mimeta$Ethnicity)
mimeta$Collection<-factor(mimeta$Collection)
design<-model.matrix(~Diagnosis+Batch+RIN+Collection+Sex+Ethnicity+AgeDeath+Hemisphere+PMI+BrainWeight+YearAutopsy,data=mimeta)
mimetacor<-cor(design[-1,-1],method='spearman')
write.csv(mimetacor,'correlation_mimeta.csv')

