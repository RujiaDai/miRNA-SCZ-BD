##differential expression by limma
library(limma)
meta$Diagnosis<-factor(meta$Diagnosis,levels=c('CTL','BP','SCZ'))
meta$RIN.squared = meta$RIN^2
meta$Batch<-factor(meta$Batch)
meta$Hemisphere<-factor(meta$Hemisphere)
meta$Sex<-factor(meta$Sex)
meta$Ethnicity<-factor(meta$Ethnicity)
meta$Collection<-factor(meta$Collection)
design<-model.matrix(~Diagnosis+Batch+RIN+Collection+Sex+Ethnicity+AgeDeath+Hemisphere+PMI+BrainWeight+YearAutopsy+sv,data=meta)
v<-voom(count.fgene.fsample,design,plot=FALSE)
fit<-eBayes(lmFit(v,design))
tt.bd<-topTable(fit,coef=2,number=Inf,sort.by='none',genelist=rownames(count.fgene.fsample))
tt.scz<-topTable(fit,coef=3,number=Inf,sort.by='none',genelist=rownames(count.fgene.fsample))
write.csv(tt.bd,'miRNA_bd_limma_dge.csv')
write.csv(tt.scz,'miRNA_scz_limma_dge.csv')

##differential expression by edgeR
library(edgeR)
exp<-DGEList(counts=count.fgene.fsample)
exp<-calcNormFactors(exp)#TMM normalization is applied to this dataset to account for compositional difference between the libraries.
# Error in glmFit.default(sely, design, offset = seloffset, dispersion = 0.05,  : 
# Design matrix not of full rank.  The following coefficients not estimable: CollectionSUN
design2<-design[,-29]
y <- estimateDisp(exp, design2, robust=TRUE)
fit <- glmQLFit(y, design2, robust=TRUE)
qlf_bd <- glmQLFTest(fit, coef=2)
qlf_scz <- glmQLFTest(fit, coef=3)
ed.scz<-topTags(qlf_scz,n=Inf)
ed.bd<-topTags(qlf_bd,n=Inf)
write.csv(ed.scz,'miRNA_scz_edger_dge.csv')
write.csv(ed.bd,'miRNA_bd_edger_dge.csv')

ed.scz<-read.csv('miRNA_scz_edger_dge.csv',header=T,sep=',',row.names=1)
ed.bd<-read.csv('miRNA_bd_edger_dge.csv',header=T,sep=',',row.names=1)

library(VennDiagram)
edger_dge<-rownames(subset(ed.scz,FDR<0.05))
limma_dge<-rownames(subset(tt.scz,adj.P.Val<0.05))
venn.diagram(
  x = list(limma_dge, edger_dge),
  category.names = c("limma" , "edgeR"),
  main='SCZ',
  filename = 'scz_dge_venn_diagramm.tiff',
  output=TRUE
)

edger_dge<-rownames(subset(ed.bd,FDR<0.05))
limma_dge<-rownames(subset(tt.bd,adj.P.Val<0.05))
venn.diagram(
  x = list(limma_dge, edger_dge),
  category.names = c("limma" , "edgeR"),
  main='BD',
  filename = 'bd_dge_venn_diagramm.tiff',
  output=TRUE
)


##volcano plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
setwd('F:/mir20200505/')
bdresult<-read.csv('miRNA_bd_edger_dge.csv',header=T,sep=',')
sczresult<-read.csv('miRNA_scz_edger_dge.csv',header=T,sep=',')

bdresult = mutate(bdresult, sig=ifelse(bdresult$FDR<0.05, "FDR<0.05", "Not sig"))
sczresult = mutate(sczresult, sig=ifelse(sczresult$FDR<0.05, "FDR<0.05", "Not sig"))
pdf("volcano_edgeR.pdf",width=600,height=800)
p1 = ggplot(bdresult, aes(logFC, -log10(FDR))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("orange", "black"))+
  labs(title = "BD")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  xlim(-2,1.5)
p2 = ggplot(sczresult, aes(logFC, -log10(FDR))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("orange", "black"))+
  labs(title = "SCZ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  xlim(-2,1.5)

grid.arrange(p1,p2,nrow=2)

#p2+geom_text_repel(data=filter(sczresult, FDR<0.05), aes(label=X))
dev.off()
 
