load('mirna289_fgene_fsample_sd3.RData')
##mrna qc
library('SummarizedExperiment')
load("mRNA/se.BrainGVEX.RData")#from PEC

datMeta<-colData(se.BrainGVEX)
v <- voom(assays(se.BrainGVEX)$counts,  plot=FALSE)
log2cpm_mrna<-v$E
mrmeta<-datMeta
save(log2cpm_mrna,mrmeta,file='mrna_qc_log2com.RData')


##integrate miRNA mRNA
rid<-paste0('X',rownames(mrmeta))

length(intersect(rid,mimeta$RNAseq.ID.1))#267
id<-intersect(rid,mimeta$RNAseq.ID.1)
mirna<-log2cpm.fgene.fsample[,match(id,mimeta$RNAseq.ID.1)]
mrna<-log2cpm_mrna[,match(id,rid)]
colnames(mirna)<-colnames(mrna)<-id
mrmeta<-mrmeta[match(id,rid),]
mimeta<-mimeta[match(id,mimeta$RNAseq.ID.1),]
save(mrna,mirna,mrmeta,mimeta,file='qc_logcpm_mrna_mirna.RData')

#load('qc_logcpm_mrna_mirna.RData')

#combat
library(sva)
library(mgcv)
library(nlme)
ibatch<-mimeta$Batch
ibatch<-ibatch[-221]
mirna<-mirna[,-221]
icombat=ComBat(dat=mirna, batch=ibatch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
rbatch<-mrmeta$final_batch[-221]
mrna<-mrna[,-221]
rcombat=ComBat(dat=mrna, batch=rbatch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
mimeta<-mimeta[-221,]
mrmeta<-mrmeta[-221,]




library(preprocessCore)
mcombine<-rbind(icombat,rcombat)
mcombine.qn<-normalize.quantiles(mcombine,copy=T)  # Quantile normalization across columns
rownames(mcombine.qn)<-rownames(mcombine)
colnames(mcombine.qn)<-colnames(mcombine)


##pca
library(ggfortify)
library(gridExtra)
pc<-prcomp(t(mcombine.qn))
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

pdf("mcombine.qn.pca.biological.pdf",width=10,height=8)
grid.arrange(g1,g3,g4,g5,g8,g10,nrow=2)
dev.off()


pdf("mcombine.qn.pca.technical.pdf",width=10,height=8)
grid.arrange(g1,g2,g6,g7,g5,g9,nrow=2)
dev.off()

##cor

#sva
library(sva)
design<-model.matrix(~Diagnosis+Collection+Sex+Ethnicity+AgeDeath+Hemisphere+BrainWeight+YearAutopsy,data=mimeta)
mod0<-mod[,-c(2:3)]
n.sv = num.sv(mcombine.qn,mod,method="leek")#0
mcom_meta<-mimeta
##variable selection
PC = prcomp(na.omit(t(scale(t(mi.mr.AllRegressed),scale=F))),scale=F) 
dM<-data.frame(mcom_meta[,c(5,7,11,12,13,14,15,16,17)])
rownames(dM)<-colnames(mcombine.qn)
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
pdf('mRNA_covariate_selection_mi.mr.AllRegressed.pdf')
ggplot(dat,aes(x=reorder(X1, value),y=value, fill=X2))+ geom_bar(stat="identity") + coord_flip() +ggtitle('Correlation with Top 5 Expression PCs') +
  xlab("") + ylab("R^2") + theme(axis.text.y=element_text(size=10), plot.title = element_text(hjust=.5))
dev.off()












mimeta$Diagnosis<-factor(mimeta$Diagnosis,levels=c('CTL','BP','SCZ'))

design<-model.matrix(~Diagnosis+Collection+PMI+Sex+Ethnicity+AgeDeath+Hemisphere+BrainWeight+YearAutopsy,data=mimeta)
Y<-mcombine.qn
X<-as.matrix(design)
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

mi.mr.AllRegressed = Y - t(X[,c(2:ncol(X))] %*% beta[c(2:nrow(beta)),])
mi.mr.keepDiagnosis = Y - t(X[,c(4:ncol(X))] %*% beta[c(4:nrow(beta)),])

save(mi.mr.AllRegressed,mi.mr.keepDiagnosis,mimeta,file='mrna_regressed.RData')
##
library(ggfortify)
library(gridExtra)
pc<-prcomp(t(mi.mr.AllRegressed))
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

pdf("mi.mr.AllRegressed.qn.pca.biological.pdf",width=10,height=8)
grid.arrange(g1,g3,g4,g5,g8,g10,nrow=2)
dev.off()


pdf("mi.mr.AllRegressed.qn.pca.technical.pdf",width=10,height=8)
grid.arrange(g1,g2,g6,g7,g5,g9,nrow=2)
dev.off()



##correlation
library(Hmisc)
mcom_reg<-mi.mr.keepDiagnosis

mcor<-rcorr(t(mcom_reg[1:666,]),t(mcom_reg[667:nrow(mcom_reg),]))

m_r<-mcor$r[667:nrow(mcom_reg),1:666]
m_p<-mcor$P[667:nrow(mcom_reg),1:666]
m_r<-melt(m_r)
m_p<-melt(m_p)
colnames(m_p)[[3]]<-'pval'
m_r$pval<-m_p$pval
nrow(subset(m_r,pval<0.05&value<0))/nrow(subset(m_r,pval<0.05))

gene<-read.csv('/zs32/home/rjdai/mirna2020/mrna_annot.csv',header=T,sep=',')
gene<-gene[,-1]
colnames(gene)[1]<-'X1'
m_r2<-merge(m_r,gene,by='X1',sort=FALSE)
m_r2$id<-paste(m_r2$X2,m_r2$name)

alltar<-read.csv('/zs32/home/rjdai/mir20200505/mRNA/all_targets.csv',header=T,sep=',')
m_r3<-m_r2[match(alltar$id,m_r2$id),]
m_r4<-m_r2[!m_r2$id%in%alltar$id,]
nrow(subset(m_r3,pval<0.05&value<0))/nrow(subset(m_r3,pval<0.05))

save(m_r2,m_r3,file='mrna_keepdiag_correlations.RData')
library(ggplot2)
library(gProfileR)
go = gprofiler(query=as.vector(unique(subset(t_r3,pval<0.05&value<0)$symbol)), max_set_size = 1000, correction_method = "fdr",hier_filtering = "strong", custom_bg = as.vector(annot$name), src_filter = c("GO"),ordered_query = F)
go = go[order(go$p.value)[1:min(10,nrow(go))],]
ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")+theme(axis.text.y = element_text(size=15))
  

pdf('mrna_target_cor_keepdiag.pdf')
plot(density(subset(m_r3,pval<0.05)$value),xlab='cor(miRNA,mRNA)',main='',ylim=c(0,8),col='red',lwd=2)
lines(density(subset(m_r4,pval<0.05)$value),col="black",add=T,lwd=2)
legend("topright",legend = c("targets","non-targets"), lty=1,lwd=2,col=c("red","black"))
dev.off()



mRNA<-subset(m_r3,pval<0.05&value<0)$id
ribo_mRNA<-subset(t_r3,pval<0.05&value<0)$id
venn.diagram(
  x = list(mRNA, ribo_mRNA),
  category.names = c("mRNA" , "ribo_mRNA"),
  main='negative correlations',
  filename = '#neg_venn_diagramm.tiff',
  cex=2,
  cat.cex=2,
  main.cex=3,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14),
  output=TRUE
)

mRNA<-subset(m_r3,pval<0.05&value>0)$id
ribo_mRNA<-subset(t_r3,pval<0.05&value>0)$id
venn.diagram(
  x = list(mRNA, ribo_mRNA),
  category.names = c("mRNA" , "ribo_mRNA"),
  main='positive correlations',
  filename = '#pos_venn_diagramm.tiff',
    cex=2,
  cat.cex=2,
  main.cex=3,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14),
  output=TRUE
)
 tidy_movies %>%
     distinct(title, year, length, .keep_all=TRUE)
mRNA_neg<-subset(m_r3,pval<0.05&value<0)$id
mRNA_pos<-subset(m_r3,pval<0.05&value>0)$id

ribo_neg<-subset(t_r3,pval<0.05&value<0)$id
ribo_pos<-subset(t_r3,pval<0.05&value>0)$id

venn.diagram(
  x = list(mRNA_neg,mRNA_pos,ribo_neg, ribo_pos),
  category.names = c("mRNA-" ,"mRNA+", "ribo-","ribo+"),
  #main='correlations',
  filename = '#sigpair_venn_diagramm.tiff',
  cex=2,
  cat.cex=1.6,
  main.cex=3,
  #cat.dist = c(0.03, 0.03,0.03,0.03),
  #cat.pos = c(-20, 14),
  output=TRUE
)

go = gprofiler(query=as.vector(m_r4[m_r4$value>0&t_r4$value<0,]$name), max_set_size = 1000, correction_method = "fdr",hier_filtering = "strong", custom_bg = as.vector(annot$name), src_filter = c("GO"),ordered_query = F)
go = gprofiler(query=as.vector(setdiff(m_r5$name,t_r5$gene)), max_set_size = 2000, correction_method = "fdr",hier_filtering = "strong", custom_bg = as.vector(annot$name), src_filter = c("GO"),ordered_query = F)
go = go[order(go$p.value)[1:min(10,nrow(go))],]
ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")+theme(axis.text.y = element_text(size=15))


mRNA<-unique(subset(m_r3,pval<0.05&value<0)$name)
ribo_mRNA<-unique(subset(t_r3,pval<0.05&value<0)$gene)
venn.diagram(
  x = list(mRNA, ribo_mRNA),
  category.names = c("mRNA" , "ribo_mRNA"),
  main='negative correlations',
  filename = '#neg_targets_venn_diagramm.tiff',
  cex=2,
  cat.cex=2,
  main.cex=3,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14),
  output=TRUE
)


mRNA<-unique(subset(m_r3,pval<0.05&value>0)$name)
ribo_mRNA<-unique(subset(t_r3,pval<0.05&value>0)$gene)
venn.diagram(
  x = list(mRNA, ribo_mRNA),
  category.names = c("mRNA" , "ribo_mRNA"),
  main='positive correlations',
  filename = '#pos_targets_venn_diagramm.tiff',
  cex=2,
  cat.cex=2,
  main.cex=3,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14),
  output=TRUE
)



mpos<-unique(subset(m_r3,pval<0.05&value>0)$name)
mneg<-unique(subset(m_r3,pval<0.05&value<0)$name)
rpos<-unique(subset(t_r3,pval<0.05&value>0)$gene)
rneg<-unique(subset(t_r3,pval<0.05&value<0)$gene)
venn.diagram(
  x = list(mpos, mneg, rneg,rpos),
  category.names = c("mRNA+","mRNA-" , "ribo-","ribo+"),
  filename = '#targets_venn_diagramm.tiff',
    cex=2,
  cat.cex=1.6,
  main.cex=3,
  #cat.dist = c(0.03, 0.03),
  #cat.pos = c(-20, 14),
  output=TRUE
)


venn.diagram(
  x = list(pos, neg),
  category.names = c("pos_targets" , "neg_targets"),
  filename = '#ribo_targets_venn_diagramm.tiff',
    cex=2,
  cat.cex=2,
  main.cex=3,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14),
  output=TRUE
)


mirnaexp1<-mcom_reg[1:666,]
sd1<-apply(mcom_reg[654:nrow(mcom_reg),],1,function(x){sd(x)/mean*(x)})
mirnaexp2<-tcom_reg[1:653,]
sd2<-apply(tcom_reg[654:nrow(tcom_reg),],1,function(x){x-mean(x)})

msdcor<-rcorr(t(mirnaexp1),sd1)
tsdcor<-rcorr(t(mirnaexp2),sd2)


q1<-t$mean1<=quantile(t$mean1)[2]
q2<-t$mean1<=quantile(t$mean1)[3]&t$mean1>quantile(t$mean1)[2]
q3<-t$mean1<=quantile(t$mean1)[4]&t$mean1>quantile(t$mean1)[3]
q4<-t$mean1>quantile(t$mean1)[4]
m1<-t[q1,]
m3<-t[q2,]
m2<-t[q2,]
m3<-t[q3,]
m4<-t[q4,]
t2<-rbind(m1,m2,m3,m4)
t2$qua<-rep(c('q1','q2','q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
t2$tar<-tmp

ggplot(t2,aes(x=qua,y=cv1,color=tar))+geom_boxplot()

t<-read.csv('ribo_mean_sd.csv',header=T,sep=',')
head(t)
annot2<-read.csv('ribo_gene.csv',header=T,sep=',')
head(annot2)
dim(annot2)
dim(t)
annot2<-read.csv('ribo_gene.csv',header=T,sep=',')
head(annot2)
dim(annot2)
dim(t)
tmp<-annot$name%in%alltar$X
tmp<-annot2$gene%in%alltar$X
table(tmp)
t.test(d$mean1[tmp],d$mean1[!tmp])
t.test(t$mean1[tmp],t$mean1[!tmp])
t.test(t$sd1[tmp],t$sd1[!tmp])
q1<-t$mean1<=quantile(t$mean1)[2]
q2<-t$mean1<=quantile(t$mean1)[3]&t$mean1>quantile(t$mean1)[2]
q3<-t$mean1<=quantile(t$mean1)[4]&t$mean1>quantile(t$mean1)[3]
q4<-t$mean1>quantile(t$mean1)[4]
m1<-t[q1,]
m3<-t[q2,]
m2<-t[q2,]
m3<-t[q3,]
m4<-t[q4,]
t2<-rbind(m1,m2,m3,m4)
head(t2)
dim(t2)
t2$name<-annot2$gene
head(t2)
t2$qua<-rep(c('q1','q2','q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
t2$tar<-tmp
head(t2)
ggplot(t2,aes(x=qua,y=cv1,color=tar))+geom_boxplot()
ggplot(t2,aes(x=qua,y=sd1,color=tar))+geom_boxplot()
t.test(subset(t2,tar=='TRUE')$sd1,subset(t2,tar!='TRUE')$sd1)
t.test(subset(t2,tar=='TRUE'&qua=='q1')$sd1,subset(t2,tar!='TRUE'&qua=='q1')$sd1)
t.test(subset(t2,tar=='TRUE'&qua=='q2')$sd1,subset(t2,tar!='TRUE'&qua=='q2')$sd1)
t.test(subset(t2,tar=='TRUE'&qua=='q3')$sd1,subset(t2,tar!='TRUE'&qua=='q3')$sd1)
t.test(subset(t2,tar=='TRUE'&qua=='q4')$sd1,subset(t2,tar!='TRUE'&qua=='q4')$sd1)
t.test(subset(t2,tar=='TRUE'&qua=='q4')$sd1,subset(t2,tar!='TRUE'&qua=='q4')$sd1)
ggplot(t2,aes(x=qua,y=sd1,color=tar))+geom_violin()
ggplot(t2,aes(x=qua,y=sd1,color=tar))+geom_boxplot()
ls()
rm(list=ls())
load('mrna_regall_correlations.RData')
d<-read.table('clipbaord',header=T,sep='\t')
d<-read.table('clipboard',header=T,sep='\t')
head(d)
dim(d)
d<-read.table('clipboard',header=T,sep='\t',row.names=1)
dim(d)
dim(d)[1]-666
dim(annot)
dim(annotalltar<-read.csv('F:/mirna2020/all_targets.csv',header=T,sep=','))
annot<-read.csv('F:/mirna2020/mrna_annot.csv',header=T,sep=',')
alltar<-read.csv('F:/mirna2020/all_targets.csv',header=T,sep=',')
annot<-read.csv('F:/mirna2020/mrna_annot.csv',header=T,sep=',')
head(annot)
dim(annot)
head(d)
q1<-d$mean1<=quantile(d$mean1)[2]
q2<-d$mean1<=quantile(d$mean1)[3]&d$mean1>quantile(d$mean1)[2]
q3<-d$mean1<=quantile(d$mean1)[4]&d$mean1>quantile(d$mean1)[3]
q4<-d$mean1>quantile(d$mean1)[4]
m1<-d[q1,]
m3<-d[q2,]
m2<-d[q2,]
m3<-d[q3,]
m4<-d[q4,]
d2<-rbind(m1,m2,m3,m4)
d2$qua<-rep(c('q1','q2'.'q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
d2$qua<-rep(c('q1','q2','q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
head(d2)
length(tmp)
d<-d[-c(1:666),]
q1<-d$mean1<=quantile(d$mean1)[2]
q2<-d$mean1<=quantile(d$mean1)[3]&d$mean1>quantile(d$mean1)[2]
q3<-d$mean1<=quantile(d$mean1)[4]&d$mean1>quantile(d$mean1)[3]
q4<-d$mean1>quantile(d$mean1)[4]
m1<-d[q1,]
m3<-d[q2,]
m2<-d[q2,]
m3<-d[q3,]
m4<-d[q4,]
d2<-rbind(m1,m2,m3,m4)
dim(d2)
head(d2)
d2$qua<-rep(c('q1','q2'.'q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
d2$qua<-rep(c('q1','q2','q3','q4'),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
head(d2)
tmp<-annot$name%in%alltar$X
length(tmp)
table(tmp)
d2$tar<-tmp
head(d2)
ggplot(t2,aes(x=qua,y=cv1,color=tar))+geom_boxplot()
ggplot(d2,aes(x=qua,y=sd1,color=tar))+geom_boxplot()
ggplot(d2,aes(x=tar,y=sd1,color=tar))+geom_boxplot()
wilcox.test(subset(d2,tar=='TRUE')$sd1,subset(d2,tar!='TRUE')$sd1)
t.test(subset(d2,tar=='TRUE')$sd1,subset(d2,tar!='TRUE')$sd1)
ggplot(d2,aes(x=tar,y=sd1,color=tar))+geom_violin()
ggplot(d2,aes(x=tar,y=sd1))+geom_violin()
library(ggpubr)
ggplot(d2,aes(x=tar,y=sd1))+geom_violin()+stat_compare_means
ggplot(d2,aes(x=tar,y=sd1))+geom_violin()+stat_compare_means()
ggplot(d2,aes(x=tar,y=sd1))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))
ggplot(d2,aes(x=tar,y=sd1))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light(base+size=15)+theme(legend.position='none')
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light(base+size=15)+theme(legend.position="none")
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light(base+size=15)
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red')
)
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light()
ggplot(d2,aes(x=tar,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light(base_size=15)+theme(legend.position='none')
ggplot(d2,aes(x=qua,y=sd1,fill=tar))+geom_violin()+stat_compare_means()+scale_x_discrete(labels=c('non-target','target'))+labs(x='',y='standard variation')+scale_fill_manual(values=c('black','red'))+theme_light(base_size=15)+theme(legend.position='none')
history(1000)


