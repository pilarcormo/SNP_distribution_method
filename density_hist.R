h<-hist(hm4$V1, breaks=10, density=10, col="lightgray", xlab="Accuracy", main="Overall") 
xfit<-seq(min(hm4$V1),max(hm4$V1),length=40) 
yfit<-dnorm(xfit,mean=mean(g),sd=sd(hm4$V1)) 
yfit <- yfit*diff(h$mids[1:2])*length(hm4$V1) 
lines(xfit, yfit, col="black", lwd=2)

