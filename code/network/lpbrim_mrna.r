load('/zs32/home/rjdai/mir20200505/mrna_adj.RData')

require(lpbrim)
N<-adj[1:664,665:nrow(adj)]
#alltar<-read.csv('/zs32/home/rjdai/mir20200505/mRNA/all_targets.csv',header=T,sep=',')
#gene1<-read.csv('/zs32/home/rjdai/mirna2020/mrna_annot.csv',header=T,sep=',')
# N<-N[,gene1$name%in%alltar$target]


N[N>=((1+.1)/2)^10]<-1
N[N<((1+.1)/2)^10]<-0
N <- N[rowSums(N)>0, colSums(N)>0] 


res = list()
for(i in seq(100)){
  res[[i]]<-findModules(N,10,sparse = FALSE)
}

a<-dim(N)[1]+dim(N)[2]
myMat = matrix(data = 0,nrow = a,ncol = a)
for(i in seq(a)){
  for(j in seq(a)){
    for(k in seq(100)){
      cls<-as.numeric(which(res[[k]]$S[i,]==1))
      myMat[i,j]<-myMat[i,j]+res[[k]]$S[j,cls]
    }
  }  
}
myMat = myMat/100
rownames(myMat)<-rownames(res[[1]]$S)
colnames(myMat)<-rownames(res[[1]]$S)

library(parallelDist)
dist.euclidean <- parDist(myMat, method = "euclidean",threads=5)
hc <- hclust(d = as.dist(dist.euclidean), method = "ave")
save(myMat,dist.euclidean,hc,file='mrna_lpbrim_100.RData')
