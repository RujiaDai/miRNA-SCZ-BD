ribo<-read.csv('ribo/14Mar19_featureCounts.csv',header=T,sep=',',row.names=1)
ribo<-ribo[,-1]
log2cpm<-voom(ribo)$E
thre<-apply(log2cpm,1,function(x){length(which(x>0.1))})>0.75*ncol(log2cpm)
log2cpm.fgene<-log2cpm[thre,]
counts.fgene<-ribo[thre,]

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
counts.fgene.fsample<-counts.fgene[,!outliers]
rmeta<-read.csv('ribo/ribo_sample2.csv',header=T,sep=',')
save(log2cpm.fgene.fsample,counts.fgene.fsample,rmeta,file='ribo_qc.RData')

#integration
load('qc_logcpm_mrna_mirna.RData')
id<-intersect(rmeta$match,colnames(mirna))#181
test1<-mirna[,match(id,colnames(mirna))]
test2<-log2cpm.fgene.fsample[,match(id,rmeta$match)]
colnames(test2)<-colnames(test1)

test1_meta<-mimeta[match(id,colnames(mirna)),]
test2_meta<-rmeta[match(id,rmeta$matched),]

library(sva)
library(mgcv)
library(nlme)
test1_meta<-test1_meta[-153,]
ibatch<-test1_meta$Batch
test1<-test1[,-153]
icombat=ComBat(dat=test1, batch=ibatch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
rbatch<-test2_meta$BatchNum[-153]
test2<-test2[,-153]
rcombat=ComBat(dat=test2, batch=rbatch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
mimeta<-test1_meta
mrmeta<-test2_meta[-153,]




library(preprocessCore)
tcombine<-rbind(icombat,rcombat)
tcombine.qn<-normalize.quantiles(tcombine,copy=T)  # Quantile normalization across columns
rownames(tcombine.qn)<-rownames(tcombine)
colnames(tcombine.qn)<-colnames(tcombine)

tcom_meta<-mimeta
tcom_meta$Diagnosis<-factor(tcom_meta$Diagnosis,levels=c('CTL','BP','SCZ'))

#sva
library(sva)
mod <-design<-model.matrix(~Diagnosis+Collection+Sex+Ethnicity+AgeDeath+Hemisphere+BrainWeight+YearAutopsy,data=tcom_meta)
mod0<-mod[,-c(2:3)]
n.sv = num.sv(tcombine.qn,mod,method="leek")#0

#tcombine.qn
PC = prcomp(na.omit(t(scale(t(mi.ri.AllRegressed),scale=F))),scale=F) 
dM<-data.frame(tcom_meta[,c(5,7,11,12,13,14,15,16,17)])
rownames(dM)<-colnames(tcombine.qn)
PC = PC$rotation[,1:5]
seqR2 = matrix(nrow=ncol(dM),ncol=5)

for(i in 1:5) {

  for(j in 1:ncol(dM)) {
    s = summary(lm(PC[,i] ~ dM[,j]))
    seqR2[j ,i] = s$adj.r.squared
  }
}

colnames(seqR2) = paste0("PC", 1:5)
rownames(seqR2) = colnames(dM)
plot(sort(rowSums(seqR2),decreasing = T))
idx=order(rowSums(seqR2),decreasing = T)
idx2 = which(rownames(seqR2)=="Diagnosis")
if(!(idx2 %in% idx)) idx = c(idx,idx2)

library(reshape); library(ggplot2)
dat = melt(seqR2[idx,])
pdf('riboRNA_covariate_selection_mi.ri.mi.ri.AllRegressed.pdf')
ggplot(dat,aes(x=reorder(X1, value),y=value, fill=X2))+ geom_bar(stat="identity") + coord_flip() +ggtitle('Correlation with Top 5 Expression PCs') +
  xlab("") + ylab("R^2") + theme(axis.text.y=element_text(size=10), plot.title = element_text(hjust=.5))
dev.off()


tcom_meta$Diagnosis<-factor(tcom_meta$Diagnosis,levels=c('CTL','BP','SCZ'))


design<-model.matrix(~Diagnosis+Collection+Sex+Ethnicity+AgeDeath+Hemisphere+BrainWeight+YearAutopsy+PMI,data=tcom_meta)
Y<-tcombine.qn
X<-as.matrix(design)
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

mi.ri.AllRegressed = Y - t(X[,c(2:ncol(X))] %*% beta[c(2:nrow(beta)),])
mi.ri.keepDiagnosis = Y - t(X[,c(4:ncol(X))] %*% beta[c(4:nrow(beta)),])
save(mi.ri.AllRegressed,mi.ri.keepDiagnosis,tcom_meta,file='ribo_regressed.RData')

##
library(ggfortify)
library(gridExtra)
pc<-prcomp(t(mi.ri.AllRegressed))
g1<-autoplot(pc, data = tcom_meta, colour = 'Diagnosis')
g2<-autoplot(pc, data = tcom_meta, colour = 'Batch')
g3<-autoplot(pc, data = tcom_meta, colour = 'AgeDeath')
g4<-autoplot(pc, data = tcom_meta, colour = 'Sex')
g5<-autoplot(pc, data = tcom_meta, colour = 'Ethnicity')
g6<-autoplot(pc, data = tcom_meta, colour = 'Collection')
g7<-autoplot(pc, data = tcom_meta, colour = 'RIN')
g8<-autoplot(pc, data = tcom_meta, colour = 'Hemisphere')
g9<-autoplot(pc, data = tcom_meta, colour = 'PMI')
g10<-autoplot(pc, data = tcom_meta, colour = 'BrainWeight')

pdf("mi.ri.AllRegressed.qn.pca.biological.pdf",width=10,height=8)
grid.arrange(g1,g3,g4,g5,g8,g10,nrow=2)
dev.off()


pdf("mi.ri.AllRegressed.qn.pca.technical.pdf",width=10,height=8)
grid.arrange(g1,g2,g6,g7,g5,g9,nrow=2)
dev.off()




library(Hmisc)
library(reshape2)
tcom_reg<-mi.ri.AllRegressed

tcor<-rcorr(t(tcom_reg[1:666,]),t(tcom_reg[667:nrow(tcom_reg),]))

t_r<-tcor$r[667:nrow(tcom_reg),1:666]
t_p<-tcor$P[667:nrow(tcom_reg),1:666]
t_r<-melt(t_r)
t_p<-melt(t_p)
colnames(t_p)[[3]]<-'pval'
t_r$pval<-t_p$pval

gene<-read.csv('/zs32/home/rjdai/mir20200505/ribo/ribo_gene.csv',header=T,sep=',')
alltar<-read.csv('/zs32/home/rjdai/mir20200505/mRNA/all_targets.csv',header=T,sep=',')


colnames(gene)[1]<-'X1'
t_r2<-merge(t_r,gene,by='X1',sort=FALSE)

t_r2$id<-paste(t_r2$X2,t_r2$gene)
t_r3<-t_r2[match(alltar$id,t_r2$id),]
t_r4<-t_r2[!t_r2$id%in%alltar$id,]

nrow(subset(t_r3,pval<0.05&value<0))/nrow(subset(t_r3,pval<0.05))
save(t_r2,t_r3,file='ribo_regall_correlations.RData')

pdf('ribo_target_regall_cor.pdf')
plot(density(subset(t_r3,pval<0.05)$value),xlab='cor(miRNA,ribo_mRNA)',main='',,ylim=c(0,8),col='red',lwd=2)
lines(density(subset(t_r4,pval<0.05)$value),col="black",add=T,lwd=2)
legend("topright",legend = c("targets","non-targets"), lty=1,lwd=2,col=c("red","black"))
dev.off()

pdf('ribo_target_regall_cor2.pdf')
plot(density(subset(t_r3,pval<0.05)$value),xlab='cor(miRNA,ribo_mRNA)',main='',,ylim=c(0,8),col='red',lwd=2)
lines(density(t_r4$value),col="black",add=T,lwd=2)
legend("topright",legend = c("targets","non-targets"), lty=1,lwd=2,col=c("red","black"))
dev.off()


pdf('ribo_target_cor.pdf')
plot(density(subset(t_r,pval<0.05)$value),xlab='cor(miRNA,ribo_mRNA)',main='',,ylim=c(0,8),col='red',lwd=2)
dev.off()
