load('mrna_regressed.RData')
##WGCNA
source('/zs32/home/rjdai/mir20200505/Hadjacency.r')
source('/zs32/home/rjdai/mir20200505/hpickSoftThreshold.r')
library(doParallel)
library(WGCNA)
options(stringsAsFactors = F)
enableWGCNAThreads(n=20)
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-hpickSoftThreshold(t(mi.mr.AllRegressed),powerVector = powers, corFnc = bicor, networkType = "hybrid2", verbose = 5)
   # Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# 1      1  0.01150  1.70          0.870 10400.0  10400.00  11500
# 2      2  0.00575 -0.57          0.866  6220.0   6250.00   7650
# 3      3  0.05770 -1.22          0.861  3790.0   3790.00   5230
# 4      4  0.14600 -1.48          0.867  2350.0   2330.00   3660
# 5      5  0.25600 -1.67          0.875  1480.0   1460.00   2630
# 6      6  0.36400 -1.74          0.892   952.0    920.00   1920
# 7      7  0.46400 -1.78          0.906   624.0    589.00   1440
# 8      8  0.53000 -1.83          0.905   418.0    382.00   1090
# 9      9  0.58200 -1.84          0.912   285.0    252.00    846
# 10    10  0.61800 -1.87          0.911   198.0    168.00    664
# 11    12  0.68600 -1.85          0.926   101.0     77.10    426
# 12    14  0.73700 -1.80          0.937    55.6     37.50    284
# 13    16  0.80900 -1.70          0.965    32.3     19.00    197
# 14    18  0.84300 -1.71          0.981    19.7     10.00    146
# 15    20  0.85300 -1.82          0.988    12.6      5.41    119




pdf(file = "mRNA_csu_powerSelect.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.8
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

adj<-Hadjacency(datExpr=t(mi.mr.AllRegressed),type='hybrid2',power=16,corFnc='bicor')
tom<-TOMsimilarity(adj, TOMType = "unsigned", verbose = 1)
save(adj,file='mrna_adj.RData')
save(tom,file='mRNA_csuwgnca_tom.RData')



##ribo
load('ribo_regressed.RData')
source('/zs32/home/rjdai/mir20200505/Hadjacency.r')
source('/zs32/home/rjdai/mir20200505/hpickSoftThreshold.r')
library(doParallel)
library(WGCNA)
options(stringsAsFactors = F)
enableWGCNAThreads(n=20)
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-hpickSoftThreshold(t(mi.ri.AllRegressed),powerVector = powers, corFnc = bicor, networkType = "hybrid2", verbose = 5)
   # Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.0683 -6.18          0.943 8030.00   8020.00 8910.0
# 2      2   0.1430 -4.43          0.942 4700.00   4670.00 5850.0
# 3      3   0.2380 -3.78          0.946 2790.00   2750.00 3930.0
# 4      4   0.3360 -3.36          0.951 1680.00   1640.00 2690.0
# 5      5   0.4460 -3.10          0.957 1020.00    987.00 1890.0
# 6      6   0.5430 -2.87          0.961  632.00    601.00 1350.0
# 7      7   0.6390 -2.78          0.967  398.00    370.00  993.0
# 8      8   0.7180 -2.77          0.972  254.00    230.00  749.0
# 9      9   0.7880 -2.73          0.979  165.00    146.00  577.0
# 10    10   0.8380 -2.71          0.981  109.00     93.70  453.0
# 11    12   0.8810 -2.65          0.972   50.50     40.40  295.0
# 12    14   0.8980 -2.59          0.962   25.00     18.20  207.0
# 13    16   0.9200 -2.44          0.962   13.30      8.47  152.0
# 14    18   0.9330 -2.28          0.962    7.49      4.12  117.0
# 15    20   0.9350 -2.14          0.958    4.48      2.07   92.4


pdf(file = "ribo_csuwgcna_powerSelect.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.8
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


adj<-Hadjacency(datExpr=t(mi.ri.AllRegressed),type='hybrid2',power=10,corFnc='bicor')
tom<-TOMsimilarity(adj, TOMType = "unsigned", verbose = 1)
save(adj,file="ribO-adj.RData")

