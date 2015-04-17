
hm<- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/hm.txt", quote="\"")
hm4<- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_4/hm.txt", quote="\"")




data <- c(hm4$V1)
length <- 20000000
h<-hist(data, breaks=50, density=50, col="gray", xlab=NULL, xlim =c(0,length), main="Homozygous SNPs on chromosome 4") 
multiplier <- h$counts / h$density
mydensity <- density(data, adjust = 2)
mydensity$y <- mydensity$y * multiplier[1]
lines(mydensity, col="darkred")


