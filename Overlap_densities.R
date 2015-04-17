
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

hm <- read.table("~/SNP_distribution_method/OF/Interesting_2/hm.txt", quote="\"")
ht <- read.table("~/SNP_distribution_method/OF/Interesting_2/ht.txt", quote="\"")


hm3 <- read.table("~/SNP_distribution_method/mob/chromosome1/hm.txt", quote="\"")
ht3 <- read.table("~/SNP_distribution_method/mob/chromosome1/ht.txt", quote="\"")


data_hm <- c(hm3$V1)
data_ht <- c(ht3$V1)
length <- 25000000

d1 <-density(data_hm, adjust = 1, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 1, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = "Chromosome 5", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = "magenta2") ##Homozygous 
lines(d2, col = "royalblue2", lty=2)

legend("topright",col=c("magenta2", "royalblue2"),lwd=1,lty=1:2,legend=c("Homozygous SNP density","Heterozygous SNP density", "k =1"), bty="n")

