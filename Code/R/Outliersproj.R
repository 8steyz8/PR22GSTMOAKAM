library(psy)
library(psych)
library(mvtnorm)
library("mvtnorm")
library("robustbase")
library("rrcov")
library(factoextra)


res.pca <- PcaHubert(microarray_imputed, scale=TRUE, k=52, kmax=52)
data_projected <- as.matrix(microarray_imputed)%*%res.pca$loadings

dim(data_projected)
p<-52

#############################################################################

##Mahalanobis distance

### Classical

Class.cov <- CovClassic(data_projected[,1:p])
plot(Class.cov)
cd <- sqrt(getDistance(Class.cov ))
c.out <- which(cd > sqrt(qchisq(0.975,p)))
length(c.out)
abline(h=sqrt(qchisq(0.975,p)),col="blue")

### MCD

MCD.cov<- CovMcd(data_projected)
plot(MCD.cov)
rd <- sqrt(getDistance(MCD.cov))
r.out <- which(rd > sqrt(qchisq(0.975,p)))
length(r.out)
abline(h=sqrt(qchisq(0.975,p)),col="blue")

#### List of the detected outliers:
cat("Robust: ", length(r.out), " outliers: ", r.out,"\n")
cat("Classical: ", length(c.out), " outliers: ", c.out,"\n")

#############################################################################

#PCA 

## Classical

pc.Class<-PcaClassic(microarray_imputed, scale=TRUE,k=52,flag=TRUE,crit=0.975)
###Distance-Distance Plot -
plot(cbind(slot(pc.Class,"sd"),slot(pc.Class,"od")),pch=19,xlab="Score distance",
     ylab="Orthogonaabline(v=slot(pc.Class,'cutoff.sd')",col="blue",lwd=2)
abline(h=slot(pc.Class,"cutoff.od"),col="blue",lwd=2)
abline(v=slot(pc.Class,"cutoff.sd"),col="blue",lwd=2)
title("Classical, k=52",cex=0.8)

cod.classical<-slot(pc.Class,"cutoff.od")
csd.classical<-slot(pc.Class,"cutoff.sd")
sd.classical<-slot(pc.Class,"sd")
od.classical<-slot(pc.Class,"od")

pca.classical.out <- which(sd.classical > csd.classical | od.classical > cod.classical)
length(pca.classical.out)


## GRID
pc.Grid<-PcaGrid(microarray_imputed, scale=FALSE,k=52,flag=TRUE,crit=0.975)
###Distance-Distance Plot -
plot(cbind(slot(pc.Grid,"sd"),slot(pc.Grid,"od")),pch=19,xlab="Score distance",
     ylab="Orthogonalabline(v=slot(pc.Grid,'cutoff.sd')",col="orange",lwd=2)
abline(v=slot(pc.Grid,"cutoff.sd"),col="orange",lwd=2)
abline(h=slot(pc.Grid,"cutoff.od"),col="orange",lwd=2)
title("PCA-Grid, k=52",cex=0.8)

cod.Grid<-slot(pc.Grid,"cutoff.od")
csd.Grid<-slot(pc.Grid,"cutoff.sd")
sd.Grid<-slot(pc.Grid,"sd")
od.Grid<-slot(pc.Grid,"od")

pca.Grid.out <- which(sd.Grid > csd.Grid| od.Grid > cod.Grid)
length(pca.Grid.out)


##ROBPCA 
pc.ROBPCA <- PcaHubert(microarray_imputed, scale=FALSE, k=52, kmax=52)
###Distance-Distance Plot -
plot(cbind(slot(pc.ROBPCA,"sd"),slot(pc.ROBPCA,"od")),pch=19,xlab="Score distance",
     ylab="Orthogoabline(v=slot(pc.ROBPCA,'cutoff.sd')",col="green3",lwd=2)
abline(v=slot(pc.ROBPCA,"cutoff.sd"),col="green3",lwd=2)
abline(h=slot(pc.ROBPCA,"cutoff.od"),col="green3",lwd=2)
title("PCA-ROBPCA, k=52",cex=0.8)

cod.rob<-slot(pc.ROBPCA,"cutoff.od")
csd.rob<-slot(pc.ROBPCA,"cutoff.sd")
sd.rob<-slot(pc.ROBPCA,"sd")
od.rob<-slot(pc.ROBPCA,"od")

pca.rob.out <- which(sd.rob > csd.rob| od.rob > cod.rob)
length(pca.rob.out)