library(psych)
microarray<-read.csv("/Users/guilhermesoares/Documents/Uni/MECD/AM/Projeto/microarray.csv")
patients<-read.csv("/Users/guilhermesoares/Documents/Uni/MECD/AM/Projeto/patients.csv")


for(i in 1:nrow(patients)) {       # for-loop over rows
  if (patients[i,'Status.at.follow.up']== "Alive") {patients[i,'Alive/Death'] <- 1}
  else{patients[i,'Alive/Death'] <- 0}
}

apply(patients[,7:11],2,summary)

boxplot(patients[,7:11])
pairs(patients[,7:11],col=patients[,13])

par(mfrow=c(2,3))
hist(patients[,7],prob=TRUE,xlab="Germinal Centre B cell signature")
hist(patients[,8],prob=TRUE,xlab="Lymph node signature")
hist(patients[,9],prob=TRUE,xlab="Proliferation signature")
hist(patients[,10],prob=TRUE,xlab="BPM6")
hist(patients[,11],prob=TRUE,xlab="MHC class II signature")

pairs.panels(patients[,7:11], smooth = FALSE, scale = FALSE, density=TRUE,
             ellipses=FALSE,digits = 2,col=patients[,13],hist.col="green")

lx<-cbind(as.vector(patients[,7:11]),as.factor(c(rep(7,nrow(patients)),rep(8,nrow(patients)),rep(9,nrow(patients)),rep(10,nrow(patients)),rep(11,nrow(patients)))))
#densityBy(lx,1,2)

par(mfrow=c(2,3))
plot(density(patients[,7]),main=colnames(patients)[7])
plot(density(patients[,8]),main=colnames(patients)[8])
plot(density(patients[,9]),main=colnames(patients)[9])
plot(density(patients[,10]),main=colnames(patients)[10])
plot(density(patients[,11]),main=colnames(patients)[11])

S=cov(patients[,7:11])
p=5
n=nrow(patients)
dif=scale(patients[,7:11],scale=FALSE)
dd=dif %*% solve(S) %*% t(dif)
d=diag(dd)
r=rank(d)
chi2q=qchisq((r-0.5)/n,p)
plot(d,chi2q,pch=20,main="",xlab="Mahalanobis distance",ylab="Chi-squared quantile",col="blue")
abline(0,1,lwd=2,col="red") #not a straight line -> not a normal distribution

library("ggplot2")
library("ggfortify")
library("gridExtra")
library("carData")
library("car")
library("factoextra")


res.pca <- prcomp(patients[,7:11], scale = TRUE)
print(res.pca)
summary(res.pca)
eig.val<-get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca)

var.pca <- get_pca_var(res.pca)
var.pca

head(var.pca$cos2)
head(var.pca$coord)

library("corrplot")
dev.off()
corrplot(var.pca$cos2, is.corr=FALSE)
fviz_cos2(res.pca, choice = "var", axes = 1:2)

fviz_pca_var(res.pca,
             col.var = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE
)

# Contributions of variables to PC1
a<-fviz_contrib(res.pca, choice = "var", axes = 1)
# Contributions of variables to PC2
b<-fviz_contrib(res.pca, choice = "var", axes = 2)
grid.arrange(a,b, ncol=2, top='Contribution of the variables to the first two PCs')

ind.pca <- get_pca_ind(res.pca)
ind.pca

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE
)

# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

autoplot(res.pca, col=(1+patients[,13]), loadings=TRUE, loadings.colour='darkorchid4', loadings.label=TRUE, loadings.label.size=3)

#proliferation signature is the best predictor of an adverse outcome


#--------------------------------ROBUST PCA---------------------------#
res.pca <- PcaHubert(microarray_imputed, scale=TRUE, k=57, kmax=62) #microarray_imputed comes from "Import Dataset" microarray_imputed.csv
data_projected <- as.matrix(microarray_imputed)%*%res.pca$loadings
print(data_projected)

write.csv(data_projected, "/Users/guilhermesoares/Desktop/data_projected.csv")
