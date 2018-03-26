## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


setwd("~/Dropbox/PaperIII/data/local_adaptation_test/LFMM/env_data")

col.rainbow <- rainbow(12)
palette(col.rainbow)


env_PC=read.table("SwAsp.Clone.Origin.data.pop.env_PC.txt",header=T)
pop=factor(env_PC$pop,levels=c("Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7","Pop8","Pop9","Pop10","Pop11","Pop12"))

png(filename="SwAsp.94samples.PC1.days5.png",width = 6, height = 6, units = 'in', res=300)
par(mar=c(4,5,1,1))
plot(env_PC$Days.plus5.degrees,env_PC$PC1,xlim=c(125,230),col=env_PC$pop,xlab=expression(paste("Number of days(>5",degree~C,")")),ylab="Environmental PC1 scores")
aggregate(env_PC[,"PC1"], list(env_PC$pop), mean)
text(137,10,"Pop11",cex=0.8)
text(155,7,"Pop9",cex=0.8)
text(155,4.5,"Pop12",cex=0.8)
text(164,2.5,"Pop10",cex=0.8)
text(175,2,"Pop7",cex=0.8)
text(184,0,"Pop8",cex=0.8)
text(194,-2.3,"Pop6",cex=0.8)
text(194,-2.8,"Pop5",cex=0.8)
text(204,-2,"Pop4",cex=0.8)
text(204,-5,"Pop3",cex=0.8)
text(214,-5.2,"Pop2",cex=0.8)
text(226,-6.8,"Pop1",cex=0.8)

z=lm(env_PC$PC1~env_PC$Days.plus5.degrees)
r=round(summary(z)$adj.r.squared,3)
#r=cor.test(env_PC$PC1,env_PC$Days.plus5.degrees)
abline(z,lwd=1,col="grey")
text(200,8,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
p=summary(z)$coefficients[,4][[2]]
ifelse(p<0.001,text(200,7,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(300,0.001,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


png(filename="SwAsp.94samples.PC2.prec.png",width = 5, height = 5, units = 'in', res=300)
plot(env_PC$Prec.Jun,env_PC$PC2,col=env_PC$pop,xlab="Prec.Jun",ylab="PC2 scores")
z=lm(env_PC$PC2~env_PC$Prec.Jun)
r=round(summary(z)$adj.r.squared,3)
#r=cor.test(env_PC$PC1,env_PC$Days.plus5.degrees)
abline(z,lwd=2)
text(45,4,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
p=summary(z)$coefficients[,4][[2]]
ifelse(p<0.001,text(45,3.3,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(300,0.001,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


plot(env_PC$Sunshine.hours.Apr.Sept,env_PC$PC3,col=env_PC$pop,xlab="Sunshine (April to Sep)",ylab="PC3 scores")


####plot of summary of PC scores on each population
env_PC1_summary=summarySE(env_PC,measurevar="PC1",groupvars="pop")
env_PC2_summary=summarySE(env_PC,measurevar="PC2",groupvars="pop")
env_PC3_summary=summarySE(env_PC,measurevar="PC3",groupvars="pop")

library(ggplot2)
ggplot(env_PC1_summary,aes(x=pop,y=PC1))+geom_errorbar(aes(ymin=PC1-se,ymax=PC1+se),width=.1)+geom_line()+geom_point()+ylab("PC1 scores")+xlab("Populations")
ggsave(filename="PC1.pop.png",height=5,width=6.5,units="in",dpi=300)
ggplot(env_PC2_summary,aes(x=pop,y=PC2))+geom_errorbar(aes(ymin=PC2-se,ymax=PC2+se),width=.1)+geom_line()+geom_point()
ggsave(filename="PC2.pop.png",height=5,width=6.5,units="in",dpi=300)
ggplot(env_PC3_summary,aes(x=pop,y=PC3))+geom_errorbar(aes(ymin=PC3-se,ymax=PC3+se),width=.1)+geom_line()+geom_point()
ggsave(filename="PC3.pop.png",height=5,width=6.5,units="in",dpi=300)

##plot for populations and days with degree higher than 5
ggplot(env_PC,aes(x=pop,y=Days.plus5.degrees))+geom_point()+geom_line()+xlab("Populations")
ggsave(filename="pop.day5.png",height=5,width=5,units="in",dpi=300)

