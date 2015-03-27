
hm1 <- read.table("~/SNP_distribution_method/OF/chromosome1/hm.txt", quote="\"")
hm2 <- read.table("~/SNP_distribution_method/OF/chromosome2/hm.txt", quote="\"")
hm3 <- read.table("~/SNP_distribution_method/OF/chromosome3/hm.txt", quote="\"")
hm4 <- read.table("~/SNP_distribution_method/OF/chromosome4/hm.txt", quote="\"")
hm5 <- read.table("~/SNP_distribution_method/OF/chromosome5/hm.txt", quote="\"")

ht1 <- read.table("~/SNP_distribution_method/OF/chromosome1/ht.txt", quote="\"")
ht2 <- read.table("~/SNP_distribution_method/OF/chromosome2/ht.txt", quote="\"")
ht3 <- read.table("~/SNP_distribution_method/OF/chromosome3/ht.txt", quote="\"")
ht4 <- read.table("~/SNP_distribution_method/OF/chromosome4/ht.txt", quote="\"")
ht5 <- read.table("~/SNP_distribution_method/OF/chromosome5/ht.txt", quote="\"")

hm <- read.table("~/SNP_distribution_method/BCF2/BCF2_chromosome3/hm.txt", quote="\"")
ht <- read.table("~/SNP_distribution_method/BCF2/BCF2_chromosome3/ht.txt", quote="\"")

shm4 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome4/hm2.txt", quote="\"")
sht4 <- read.table("~/SNP_distribution_method/Sch/samse_chromosome4/ht2.txt", quote="\"")

hm_copy <- read.table("~/SNP_distribution_method/OF/chromosome2/hm_copy.txt", quote="\"")
hm_copy <- read.table("~/SNP_distribution_method/OF/chromosome2/ht_copy.txt", quote="\"")
data_hm <- c(hm_copy$V1)
data_ht <- c(ht_copy$V1)
length <- 20000000

d1 <-density(data_hm, adjust = 1, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 1, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = "Chromosome 2", xlim =c(10000000,length), xlab = "Chromosome length", ylab ="Frequency")
lines(d1, col = "magenta2") ##Homozygous 
lines(d2, col = "royalblue2", lty=2)

legend("topright",col=c("magenta2", "royalblue2"),lwd=1,lty=1:3,legend=c("Homozygous SNP density","Heterozygous SNP density", "k =1"), bty="n")

