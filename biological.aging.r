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
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
big.ba.measures[,1]<-Table[,1]
names<-colnames(Table)
colnames(big.ba.measures)<-c(names[1],"MLR")
table <- Table[,3:ncol(Table)]
age <- table$age
fit <- lm(age ~ ., data=table) 
BA_mlr0 <- fitted(fit)
BA_mlr <- scale(BA_mlr0, center=TRUE, scale=TRUE)*sd(Table[,3]) + mean(Table[,3])
big.ba.measures[,2]<-BA_mlr
return(big.ba.measures)
}

###BA by JW
BA.by.JW<-function(Table){
ba.measures<-data.frame(matrix(ncol=length(colnames(Table)),nrow=length(Table[,1])))
ba.measures.variances<-data.frame(matrix(ncol=length(colnames(Table)),nrow=length(Table[,1])))
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
for (i in 1:3){
	ba.measures[,i]<-Table[,i]
	ba.measures.variances[,i]<-Table[,i]
}
big.ba.measures[,1]<-Table[,1]
colnames(ba.measures)<-colnames(Table)
colnames(ba.measures.variances)<-colnames(Table)
names<-colnames(Table)
colnames(big.ba.measures)<-c(names[1],"JW")
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
big.ba.measures[i,2]<-BAJW
}
holder<-big.ba.measures[,2]
holder2<-scale(holder, center=TRUE, scale=TRUE)*sd(Table[,3]) + mean(Table[,3])
big.ba.measures[,2]<-holder2
return(big.ba.measures)
}

###BA by TLS or (multivariate) orthogonal regression
BA.by.TLS<-function(Table){
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
big.ba.measures[,1]<-Table[,1]
names<-colnames(Table)
colnames(big.ba.measures)<-c(names[1],"TLS")
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
big.ba.measures[,2]<-BA_tls
return(big.ba.measures)
}

###Function to calculate biological age using Bayesian Model Averaging
BA.by.BMA<-function(Table){
library(BMA)
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
big.ba.measures[,1]<-Table[,1]
names<-colnames(Table)
colnames(big.ba.measures)<-c(names[1],"BMA")
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
big.ba.measures[,2]<-BA_BMA
return(big.ba.measures)
}

###Function to calculate biological age using the continuous frailty index
BA.by.CFI<-function(Table){
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
big.ba.measures[,1]<-Table[,1]
names<-colnames(Table)
colnames(big.ba.measures)<-c(names[1],"CFI")
pbig <- NCOL(Table)
age <- Table[,3]
met_syn_mat1 <- Table[,4:pbig]
svd_met_syn <- svd(met_syn_mat1)
u1 <- svd_met_syn$u 
corrv <- cor(u1,age)
ind_v <- as.vector(ifelse(corrv<0,-1,1))
ind_mat <- diag(ind_v)
u2 <- svd_met_syn$u%*%ind_mat
cFI <- rowSums(u2)
fit2 <- lm(cFI~age)
summary2 <- summary(fit2)
alpha <- summary2$coefficients[1,1]
beta <- summary2$coefficients[2,1]
BA_cfi0 <- (cFI-alpha)/beta
BA_cfi <- scale(BA_cfi0, center=TRUE, scale=TRUE)*sd(age)+ mean(age)
big.ba.measures[,2]<-BA_cfi
return(big.ba.measures)
}

###Function to calculate biological age using the KD method
BA.by.KD<-function(TABLE){
big.ba.measures<-data.frame(matrix(ncol=2,nrow=length(TABLE[,1])))
big.ba.measures[,1]<-TABLE[,1]
names<-colnames(TABLE)
colnames(big.ba.measures)<-c(names[1],"KD")
names_TABLE<- as.vector(names(TABLE))
AGE<- TABLE$age
i=4
BAKD1<-rep(0,dim(TABLE)[1])
demo_vec<-rep(0,dim(TABLE)[1])
while(i<= dim(TABLE)[2]){
	bio_marker<- names(TABLE)[i]
	BIO <- TABLE[,i]
	REGRES<- lm(AGE~BIO)
	RESULTS<- summary(REGRES)
	intercept<- RESULTS$coefficients[1,1]
	slope<- RESULTS$coefficients[2,1]
	MARKER<- names(BIO)
	RMSE<- RESULTS$sigma
	s2<- RMSE^2
	BAKD1_Single <- as.vector(   ((BIO - intercept)*(slope/s2)))
	dem<- (slope/sqrt(s2))^2 
	demo_vec<- demo_vec + dem 
	BAKD1<- BAKD1 + BAKD1_Single		
	i=i+1
}
BAKD1<- BAKD1/demo_vec
sum_sq_r2<- 0
sum_r2 <- 0 
i=4
while(i<= dim(TABLE)[2]){
	bio_marker<- names(TABLE)[i]
	BIO <- TABLE[,i]
	REGRES<- lm(AGE~BIO)
	RESULTS<- summary(REGRES)
	r2<- RESULTS$r.squared		
	sum_sq_r2<- sum_sq_r2 + (r2/sqrt(1 - r2))
	sum_r2<- sum_r2 + (sqrt(r2) /sqrt(1 - r2))
	i=i+1
}
rchar<- sum_sq_r2/sum_r2
s_diff<- (sqrt(1-(rchar)^2)/rchar) * ((max(AGE) - min(AGE))/sqrt(12*(dim(TABLE)[2])))
s2b<-   (sum(((BAKD1 - AGE) -  (1/dim(TABLE)[1]) * sum(BAKD1 - AGE))^2)/ dim(TABLE)[1]) - s_diff^2
BAKD2<-rep(0,dim(TABLE)[1])
i=4
demo_vec<-rep(0,dim(TABLE)[1])
while(i<= dim(TABLE)[2]){
	bio_marker<- names(TABLE)[i]
	BIO <- TABLE[,i]
	REGRES<- lm(AGE~BIO)
	RESULTS<- summary(REGRES)
	intercept<- RESULTS$coefficients[1,1]
	slope<- RESULTS$coefficients[2,1]
	MARKER<- names(BIO)
	AGE_IND<- AGE[i]
	RMSE<- RESULTS$sigma
	s2<- RMSE^2
	BAKD2_Single <- as.vector(     ((BIO - intercept)*(slope/s2))     + AGE_IND/s2b )
	demo_vec<-demo_vec +      ((slope/sqrt(s2))^2 + 1/s2b   )
	BAKD2<- BAKD2 + BAKD2_Single		
	i=i+1
}
BAKD2<- BAKD2/demo_vec
BA_kd0 <- BAKD2
BA_kd <- scale(BA_kd0, center=TRUE, scale=TRUE)*sd(AGE) + mean(AGE)
big.ba.measures[,2]<-BA_kd
return(big.ba.measures)
}

calculate.ba.measures<-function(Table){
Newtable<-data.frame(matrix(ncol=2,nrow=length(Table[,1])))
Newtable[,1]<-Table[,1]
Newtable[,2]<-Table[,3]
names<-colnames(Table)
colnames(Newtable)<-c(names[1],"CA")

MLR<-BA.by.MLR(Table)
JW<-BA.by.JW(Table)
TLS<-BA.by.TLS(Table)
KD<-BA.by.KD(Table)
CFI<-BA.by.CFI(Table)
BMA<-BA.by.BMA(Table)
NEW0<-merge(Newtable,MLR,by=names[1])
NEW1<-merge(NEW0,JW,by=names[1])
NEW2<-merge(NEW1,TLS,by=names[1])
NEW3<-merge(NEW2,KD,by=names[1])
NEW4<-merge(NEW3,CFI,by=names[1])
NEW5<-merge(NEW4,BMA,by=names[1])
write.table(NEW5, file="BA_table_PC_scores_met_syn1.csv", sep=",",row.names=FALSE)
return(NEW5)
}