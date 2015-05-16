


hm1 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/hm.txt", quote="\"")
ht1 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/ht.txt", quote="\"")
hm4 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup2_2_4/hm.txt", quote="\"")
ht4 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup2_2_4/ht.txt", quote="\"")



data_hm <- c(hm1$V1)      
data_ht <- c(ht1$V1)
length <- 31000000
 
d1 <-density(data_hm, adjust = 1, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 1, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = "Chromosome 1", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = "magenta2") ##Homozygous 
lines(d2, col = "royalblue2", lty=2)

legend("topleft",col=c("magenta2", "royalblue2"),lwd=1,lty=1:2,
       legend=c("Homozygous SNP density","Heterozygous SNP density", "k =1"), bty="n")




