p_hm1 <- read.table("~/SNP_distribution_method/Sch/chromosome1/hm.txt", quote="\"")


hm1 <- read.table("~/SNP_distribution_method/OF/chromosome1/hm.txt", quote="\"")
hm2 <- read.table("~/SNP_distribution_method/OF/chromosome2/hm.txt", quote="\"")
hm3 <- read.table("~/SNP_distribution_method/OF/chromosome3/hm.txt", quote="\"")
hm4 <- read.table("~/SNP_distribution_method/OF/chromosome4/hm.txt", quote="\"")
hm5 <- read.table("~/SNP_distribution_method/OF/chromosome5/hm.txt", quote="\"")

shm1 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome1/hm2.txt", quote="\"")
shm2 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome2/hm2.txt", quote="\"")
shm3 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome3/hm2.txt", quote="\"")
shm4 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome4/hm2.txt", quote="\"")
shm5 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome5/hm2.txt", quote="\"")

hm <- read.table("~/SNP_distribution_method/BCF2/chromosome3/hm2.txt", quote="\"")
ht <- read.table("~/SNP_distribution_method/BCF2/chromosome3/ht2.txt", quote="\"")

data <- c(shm4$V1)
length <- 25000000
h<-hist(data, breaks=50, density=50, col="gray", xlab=NULL, xlim =c(0,length), main="Homozygous SNPs on chromosome 5") 
multiplier <- h$counts / h$density
mydensity <- density(data, adjust = 2)
mydensity$y <- mydensity$y * multiplier[1]
lines(mydensity, col="darkred")


