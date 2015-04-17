
hm2 <- read.table("~/SNP_distribution_method/BCF2/BCF2_chromosome3/hm.txt", quote="\"")
ht2 <- read.table("~/SNP_distribution_method/BCF2/BCF2_chromosome3/ht.txt", quote="\"")

hm3 <- read.table("~/SNP_distribution_method/BCF2_parent/chromosome3/hm.txt", quote="\"")
ht3 <- read.table("~/SNP_distribution_method/BCF2_parent/chromosome3/ht.txt", quote="\"")

hm <- read.table("~/SNP_distribution_method/OF/Interesting_2/hm.txt", quote="\"")
data <- c(hm$V1)
length <- 25000000
h<-hist(data, breaks=50, density=50, col="gray", xlab=NULL, xlim =c(10000000,length), main="Homozygous SNPs on chromosome 2") 
multiplier <- h$counts / h$density
mydensity <- density(data, adjust = 2)
mydensity$y <- mydensity$y * multiplier[1]
lines(mydensity, col="darkred")


