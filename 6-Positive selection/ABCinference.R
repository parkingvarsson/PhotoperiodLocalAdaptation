#######################################################################
# Joint inference of s and T assuming an equilibrium population
# Example script for mouse coat color, modified from scrip provided in Ormond et al. (2016).

#There are four sections in this script
# 1) Run the simulations for ABC inference
# A file containing 500,000 simulations for this example can be read in at the end of part 1
# 2) Calculate PLS components from the statistics
# 3) Run simulations of 100 pseudo observables to establish how well the method works across different orders of magnitude for s and T
# 4) Apply the method to the aspen data


#install and open the required R packages
library("pls"); library ("abc"); library("MASS");library("data.table");library("RColorBrewer")

#1) Run the simulations for ABC inference
#Identify the parameters for input to the simulations to be run in msms (please refer to msms manual)
#here all parameters are from Linnen et al. 2013.  A sequence length of 40kb either side of the serine deletion was used. After filtering for genotype quality a sample size of 48 was obtained.
#the length of the sequence
L=25000
#the effective population size
Ne=92080
#the sample size
samplesize=52
#the number of replicates (here we use one replicate per draw of s and T)
nrep=1
#the position of the selected mutation halfway along the sequence
Sp=0.5
#Assigne values to theta and rho
#theta
mu=3.75*10^-8
theta=4*Ne*mu*L
#rho
r= 0.729*10^-8
rho=4*Ne*r*L


#set the filename
run=1  # different run numbers can be used to generate parallel sets of simulations
filename <- paste('SwAsp',run,sep='')

#run the simulations
nsim <-200000 # set the number of simulations.  There is a file with 500,000 simulations to be read in at the end of part 1.
msstats=c()
for (i in 1:nsim)  {
  a <- runif(1, -4, -0.5) # draw the selection coefficient s from a log uniform prior log10(T)~U(-4;-0.5)
  selec_coef=10^a
  alpha_homo=2*Ne*selec_coef  #set alpha for the homozygote and heterozygote genotypes using the selection coefficient
  alpha_hetero=(alpha_homo)/2
  b <- runif(1,-4,-0.5)  #draw T from from a log uniform prior log10(T)~U(-4;-0.5)
  T<-10^b
  # write these parameters to a file for input to the msms command line and run the msms simulation
  write.table(data.frame("alpha_homo"=format(alpha_homo,scientific=F),"alpha_hetero"=format(alpha_hetero,scientific=F),"T"=format(T,scientific=F)), file=paste('Inputs_', filename,'.txt',sep=''), quote=F,row.name=F, col.name=F)
  String <-system(paste("~/src/msms/bin/msms 52 1 -N 92000 -t 345.3 -r 67.12632 -SAA tbs -SAa tbs -Sp 0.5 -SF tbs -SFC < Inputs_",filename,".txt",sep=''), intern=TRUE)
  #write the SNP output from the msms simulation to a file
  write (String, file=paste("string_",filename,".txt",sep=''), sep ="\t")
  #calculate statistics from the msms output using mssstats and write these to a file
  system (paste("cat string_",filename,".txt | msstats >  msstats_",filename,".txt",sep=''))
  stats <-read.table(file=paste("msstats_",filename,".txt",sep=''), header = TRUE)
  #combine the calculated statistics with the values of s and T and record these in the object msstatstotal
  stats2 <- cbind(selec_coef,T,stats)
  msstats <- rbind(msstats,stats2)
}
write.table (msstats, file=paste("msstats_SwAsp",filename,".txt",sep=''), sep ="\t")

#Here we read in a file from previously run simulations
msstats<-fread(file="equilibrium_mouse_sims.txt", header=TRUE)

#2) Calculate PLS components from the statistics

#check that all statistics have been calculated correctly
row.has.na <- apply(msstats, 1, function(x){any(is.na(x))})
sum(row.has.na) #ideally this should be 0 so that the log priors are uniformly sampled.
msstats <- msstats[!row.has.na,]

#calculate PLS loadings based on a subset of the first 10,000 simulations (approach follows Wegmann et al. 2009)
a<- msstats[1:10000,]
stat<-a[,3:23]; param<-a[,1:2]; #the number of segregating sites S is excluded from the calculation as we have conditioned on S in the msms simulations
#remove invariant statistics
stat<-stat[,-14]
##stat<-stat[,-3]

#standardize the parameters and the statistics
mymeanparam <- c()
mysdparam <- c()
for(i in 1:length(param)){
  mymeanparam <- c(mymeanparam, mean(param[,i])); mysdparam <- c(mysdparam, sd(param[,i]));
  param[,i]<-(param[,i]-mean(param[,i]))/sd(param[,i]);}
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stat)){
  myMax<-c(myMax, max(stat[,i])); myMin<-c(myMin, min(stat[,i]));
  stat[,i]<-1+(stat[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
# apply Box-Cox transformation to normalize the statistics prior to PLS calculation
for(i in 1:length(stat)){
  d<-cbind(stat[,i], param);
  mylm<-lm(as.formula(d), data=d);
  myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
  lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
  myGM<-c(myGM, mean(exp(log(stat[,i]))));
}
#standardize the Box-Coxed statistics
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stat)){
  stat[,i]<-(stat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  myBCSDs<-c(myBCSDs, sd(stat[,i]));
  myBCMeans<-c(myBCMeans, mean(stat[,i]));
  stat[,i]<-(stat[,i] -myBCMeans[i])/myBCSDs[i];
}
#Calculate the PLS components
myPlsr<-plsr(as.matrix(param)~as.matrix(stat), scale=F, validation='LOO');
#write the PLS information to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:10){myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i])}
write.table(cbind(names(stat), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs,myPlsrDataFrame), file="plssummary_SwAsp.txt", sep="\t", col.names=F,row.names=F, quote=F);
# the reduction in RMSE for each parameter from using PLS can be visualised
plot(RMSEP(myPlsr));
# the component loadings can be generated using
#loadings(myPlsr)

#the required PLS information for the next step - the ABC inference - can be recovered from the summary file as follows:
#summary <- read.table("plssummary_mice.txt", header=FALSE)
#myMax <- summary[,2]
#myMin <- summary[,3]
#lambda <- summary[,4]
#myGM <- summary[,5]
#myBCMeans <- summary[,6]
#myBCSDs <- summary[,7]
#myPlsrDataFrame <- summary[,8:length(summary)]


#convert the statistics from the simulations into PLS components
b<-msstats
scores<-b[,3:23]; param<-b[,1:2];
#remove invariant statistics (same steps as for PLS calculation)
scores <- scores[,-14]
#scores <- scores[,-3]
#standardize the statistics
for(i in 1:length(scores)){
  scores[,i]<- 1+(scores[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
#apply the Box Cox transformation and normalize
for(i in 1:length(scores)){
  scores[,i]<-(scores[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  scores[,i]<-(scores[,i] -myBCMeans[i])/myBCSDs[i];
}
#convert into PLS components (here the top 7 components are used)
scores <- (as.matrix(scores))%*%(as.matrix(myPlsrDataFrame[,1:7]))

#3) Run simulations of 100 pseudo observables to establish how well the method works across different orders of magnitude for s and T
#E.g. true values of s=0.1, T=0.01
s_true=0.1
alpha_homo=2*s_true*Ne
alpha_homo
alpha_hetero=alpha_homo/2
alpha_hetero
T_true=0.01

#redefine the mode function
mode=function(v) {
  dens=density(v)
  mode=dens$x[which.max(dens$y)]
}

#run the inference for the pseudo-observables and calculate relative bias and RMSE
#input the values of s_true and T_true into the msms command line
ztot=matrix(0,nrow=100,ncol=100)
s_box=c()
T_box=c()
for (i in 1:100) {
  #run the pseudo observable simulation and convert calculated statistics into scores applying the same steps as above
  system("msms 48 1 -N Ne -t theta -r rho -s 842 -SAA 10616 -SAa 5308 -Sp 0.5 -SF 0.01 -SFC | msstats >test.txt")
  test<-read.table("test.txt", header = TRUE)
  test <- test[,-1]
  test <- test[,-13]
  test <- test[,-3]
  for(i in 1:length(test)){
    test[,i]<- 1+(test[,i] -myMin[i])/(myMax[i]-myMin[ i]);
  }
  for(i in 1:length(test)){
    test[,i]<-(test[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
    test[,i]<-(test[,i] -myBCMeans[i])/myBCSDs[i];
  }
  newtest <- (as.matrix(test))%*%(as.matrix(myPlsrDataFrame[,1:7]))
  #run the ABC inference using a tolerance of 0.005 and rejection ABC
  sim<-abc(target = c(newtest[1,c(1,2,3,4,5,6,7)]), param=param[,c(1,2)], sumstat=scores[,c(1,2,3,4,5,6,7)], method ="rejection", tol=0.005, transf="none")
  post_s=(sim$unadj.values[,1])
  post_T=(sim$unadj.values[,2])
  #calculate the joint posterior density and record it for the cumulative joint density plots
  z <- kde2d(log10(post_s),log10(post_T),lims=c(-4,0,-4,0),n=100)
  ztot=ztot+z$z
  maxindex=which(z$z==max(z$z),arr.ind=T)
  #calculate point estimates of s and T for each pseudo-observable simulation and record these
  s_pred <- z$x[maxindex[1]]
  T_pred <- z$y[maxindex[2]]
  s_box=c(s_box,s_pred)
  T_box=c(T_box,T_pred)
}

#Display the cumulative joint density plots from the simulations
par(mfrow=c(1,1))
image(z$x,z$y,ztot,xlab="log(s)",ylab="log(T)", main="Cumulative joint distribution", cex.main=1, cex.lab=1)
points(x=log10(s_true),y=log10(T_true),pch=4, cex=2)

#Calculate relative bias and RMSE
RMSE_s=sqrt((sum((10^s_box)-s_true)^2)/100)
RMSE_T=sqrt((sum((10^T_box)-T_true)^2)/100)
relbias_s=(mean((10^s_box)-s_true))/s_true
relbias_T=(mean((10^T_box)-T_true))/T_true
relbias_s
RMSE_s
relbias_T
RMSE_T

#4) Apply the method to the mouse data

#calculate statistics from the SNP mouse data
##system (("cat mousedata80kb.txt | msstats > msstatstest.txt"), intern=TRUE)

test <-read.table("chr10.msstats.txt", header = TRUE)
#remove S and invariant statistics and convert the statistics using the PLS loadings (same steps as above)
test <- test[,-c(14)]
#test <- test[,-14]
#test <- test[,-3]
for(i in 1:length(test)){
  test[,i]<- 1+(test[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
for(i in 1:length(test)){
  test[,i]<-(test[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
  test[,i]<-(test[,i] -myBCMeans[i])/myBCSDs[i];
}
newtest <- (as.matrix(test))%*%(as.matrix(myPlsrDataFrame[,1:7]))

# Do the ABC inference
sim<-abc(target = c(newtest[1,c(1,2,3,4,5,6,7)]), param=param, sumstat=scores[,c(1,2,3,4,5,6,7)], method ="rejection", tol=0.005, transf="none")
post_s=(sim$unadj.values[,1])
post_T=(sim$unadj.values[,2])
z <- kde2d(log10(post_s),log10(post_T),lims=c(-4,0,-4,0),n=100)
maxindex=which(z$z==max(z$z),arr.ind=T)
s_pred <- z$x[maxindex[1]]
T_pred <- z$y[maxindex[2]]

#Display the joint posterior density plot with a point estimate from the mode (indicated by a black cross)
pdf("TMRCA_s_northern_allele.pdf",width=5,height=5)
colfunc<-colorRampPalette(brewer.pal(9,"OrRd"))
par(mfrow=c(1,1))
par(las=1)
image(x=z$x,y=z$y,z=z$z,xlab="log(s)",ylab="log(T)",cex.main=1.3, cex.lab=1.3, main="",col = colfunc(100))
points(x=s_pred,y=T_pred, pch=4, cex=1.2)
dev.off()
#Point estimates for s and T
Est_s<-10^(s_pred)
Est_T<-10^(T_pred)
Est_s
Est_T

#Generate 95% confidence intervals
quantile(post_s,probs=c(0.025,0.975))
quantile(post_T,probs=c(0.025,0.975))

#Display histograms of the marginal posterior density
par(mfrow=c(1,2))
hist(log10(post_s), xlim=c(-4,0), main="Predicted s", xlab="log(s)", breaks=20)
abline(v=s_pred,b=0, col="red", cex=100)
hist(log10(post_T),xlim=c(-4,0),main="Predicted T", xlab="log(T)", breaks=20)
abline(v=T_pred,b=0, col="red", cex=100)
