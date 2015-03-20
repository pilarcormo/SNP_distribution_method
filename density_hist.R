data <- c(hm5$V1)
length <- 25000000
h<-hist(data, breaks=50, density=50, col="gray", xlab="Chromosome length", xlim =c(0,length), main="Homozygous SNPs on chromosome 5") 
multiplier <- h$counts / h$density
mydensity <- density(data, adjust = 2)
mydensity$y <- mydensity$y * multiplier[1]
lines(mydensity, col="darkred")

