###Function to regress out sex
regress.sex<-function(Table){
regress.out.sex<-data.frame(matrix(ncol=length(colnames(Table)),nrow=length(Table[,1])))
for (i in 1:3){
	regress.out.sex[,i]<-Table[,i]
}
colnames(regress.out.sex)<-colnames(Table)
for( i in 4:length(colnames(Table))){
	Result<-NULL
	Result<-lm(Table[,i]~Table[,2],na.action=na.exclude)
	regress.out.sex[names(Result$residuals),i]<-Result$residuals
}
return(regress.out.sex)
}

###Function to perform principal components analysis using the singular value decomposition for a given dataset
calculate.pca.using.svd<-function(table,cut_off){
pca.using.svd<-data.frame(matrix(ncol=3,nrow=length(table[,1])))
for (i in 1:3){
	pca.using.svd[,i]<-table[,i]
}
colnames(pca.using.svd)<-colnames(table[,1:3])
p0 <- NCOL(table)
ms <- table[,4:p0]
svdms <- svd(ms)
u <- svdms$u
d <- svdms$d
v <- svdms$v
Eigenvalues <- d^2
sum_Eigenvalues <- sum(Eigenvalues)
Variance <- Eigenvalues/sum_Eigenvalues
Cumulative_variance <- cumsum(Variance)
eigvar <- rbind(Eigenvalues,Variance,Cumulative_variance)
ms1 <- Cumulative_variance >= cut_off
ms2 <- Cumulative_variance[ms1]
ms3 <- ms2[1]
p <- which(Cumulative_variance==ms3)
PC_scores <- u[,1:p]
report <- c(p, ms3)
plot(seq(1:length(Eigenvalues)), Eigenvalues, xlab="Eigenvalue order", ylab="Eigenvalue",main="Scree plot", typ="l",col="blue",lwd=3)
write.table(eigvar, file="Eig_var_Cumulative_variance.csv", sep=",")
PC_scores2 <- merge(pca.using.svd, PC_scores, by = "row.names", all = TRUE)
write.table(PC_scores, file="PC_scores.csv", sep=",")
PC_scores3 <- PC_scores2[,2:ncol(PC_scores2)]
return(PC_scores3)
}

#########
#########
#BA measure estimation approaches
#########
#########

###Classical biological age (BA) by multiple linear regression
BA.by.MLR<-function(Table){
big.ba.measures<-data.frame(matrix(ncol=4,nrow=length(Table[,1])))
for (i in 1:3){
	big.ba.measures[,i]<-Table[,i]
}
colnames(big.ba.measures)<-c(colnames(Table[,1:3]),"MLR")
table <- Table[,3:ncol(Table)]
age <- table$age
fit <- lm(age ~ ., data=table) 
BA_mlr0 <- fitted(fit)
BA_mlr <- scale(BA_mlr0, center=TRUE, scale=TRUE)*sd(Table[,3]) + mean(Table[,3])
big.ba.measures[,4]<-BA_mlr
return(big.ba.measures)
}

###BA by JW
BA.by.JW<-function(Table){
ba.measures<-data.frame(matrix(ncol=length(colnames(Table)),nrow=length(Table[,1])))
ba.measures.variances<-data.frame(matrix(ncol=length(colnames(Table)),nrow=length(Table[,1])))
big.ba.measures<-data.frame(matrix(ncol=4,nrow=length(Table[,1])))
print(dim(ba.measures))
print(dim(ba.measures.variances))
print(dim(big.ba.measures))
flush.console()
for (i in 1:3){
	ba.measures[,i]<-Table[,i]
	ba.measures.variances[,i]<-Table[,i]
	big.ba.measures[,i]<-Table[,i]
}
colnames(ba.measures)<-colnames(Table)
colnames(ba.measures.variances)<-colnames(Table)
colnames(big.ba.measures)<-c(colnames(Table[,1:3]),"BAJW")
for( i in 4:length(colnames(Table))){
	Result<-NULL
	Result_summary<-NULL
	alpha<-NULL
	beta<-NULL
	residual_sigma<-NULL
	V_ba_i<-NULL
	Result<-lm(Table[,i]~Table[,3])
	Result_summary<-summary.lm(Result)
	alpha<-Result_summary$coefficients[1,1]
	beta<-Result_summary$coefficients[2,1]
	residual_sigma<-Result_summary$sigma
	ba.measures[,i]<-(Table[,i]-alpha)/beta
	ba.measures.mean<-mean(ba.measures[,i],na.rm=TRUE)
	ba_var<-((residual_sigma^2)/(beta^2))*(1+(1/length(Result$residuals))+(((ba.measures[,i] - ba.measures.mean)^2)/((sum((Table[names(Result$residuals),3])^2))-((1/length(Table[names(Result$residuals),3]))*(sum(Table[names(Result$residuals),3])^2)))))
	ba.measures.variances[,i]<-ba_var
}
for( i in 1:length(Table[,1])){
BAJW<-(sum(ba.measures[i,4:length(colnames(Table))]/ba.measures.variances[i,4:length(colnames(Table))])/sum(1/ba.measures.variances[i,4:length(colnames(Table))]))
big.ba.measures[i,4]<-BAJW
}
holder<-big.ba.measures[,4]
holder2<-scale(holder, center=TRUE, scale=TRUE)*sd(Table[,3]) + mean(Table[,3])
big.ba.measures[,4]<-holder2
return(big.ba.measures)
}

###BA by TLS or (multivariate) orthogonal #regression
BA.by.TLS<-function(Table){
big.ba.measures<-data.frame(matrix(ncol=4,nrow=length(Table[,1])))
for (i in 1:3){
	big.ba.measures[,i]<-Table[,i]
}
colnames(big.ba.measures)<-c(colnames(Table[,1:3]),"TLS")
pbig <- NCOL(Table)
met_syn_mat0 <- Table[,3:pbig]
met_syn_mat1 <- as.matrix(met_syn_mat0)
age <- met_syn_mat1[,1] 
pmat1 <- NCOL(met_syn_mat1)
met_syn_mat <- cbind(met_syn_mat1[,2:pmat1],age)
met_syn_mat2 <- met_syn_mat1[,2:pmat1]
svd_met_syn <- svd(met_syn_mat)
v <- svd_met_syn$v
n <- NROW(v)
p <- NCOL(v)
v11 <- v[1:(n-1),1:(p-1)]
v12 <- v[1:(n-1),p]
v21 <- v[n,1:(p-1)]
v22 <- v[n,p]
beta_v <- -v12%*%solve(v22)
BA_tls0 <- met_syn_mat2%*%beta_v
BA_tls <- scale(BA_tls0, center=TRUE, scale=TRUE)*sd(age) + mean(age)
big.ba.measures[,4]<-BA_tls
return(big.ba.measures)
}

###Function to calculate biological age using Bayesian Model Averaging
calculate.ba.measures.BMA<-function(Table){
library(BMA)
big.ba.measures<-data.frame(matrix(ncol=4,nrow=length(Table[,1])))
for (i in 1:3){
	big.ba.measures[,i]<-Table[,i]
}
colnames(big.ba.measures)<-c(colnames(Table[,1:3]),"BMA")
p <- NCOL(Table)
ms1 <- Table[,4:p]
age <- Table[,3]
msbicreg <- bicreg(ms1,age)
summ <- summary(msbicreg)
ones <- rep(1,times=NROW(ms1))
ms2 <- as.matrix(cbind(ones,ms1))
beta_v <- as.numeric(summ[1:(p-2),2])
predms <- ms2%*%beta_v
BA_BMA <- scale(predms, center=TRUE, scale=TRUE)*sd(age) + mean(age)
big.ba.measures[,4]<-BA_BMA
return(big.ba.measures)
}





###Example
##BioAgeData <- read.csv(file="phensForScreen.csv",head=T)
##Sex.corrected<-regress.sex(BioAgeData)
##PCA.of.data<-calculate.pca.using.svd
##ba.measures.measured<-calculate.ba.measures(PCA.of.data)
##ba.via.bma<-calculate.ba.measures.BMA(PCA.of.data)

