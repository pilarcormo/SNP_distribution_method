


hm1 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/hm.txt", quote="\"")
ht1 <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/ht.txt", quote="\"")
hm4 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/hm_nocen.txt", quote="\"")
ht4 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/ht_nocen.txt", quote="\"")




hm <- read.table("~/SNP_distribution_method/Reads/Alpina/arabis_bc1f2_1/Interesting_bc1f2/hm2.txt", quote="\"")
ht <- read.table("~/SNP_distribution_method/Reads/Alpina/arabis_bc1f2_1/hm.txt", quote="\"")
#ratio <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/ratios_BC.txt", quote="\"")
ratio <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/sup1_nocen_chr4_100kb/SDM_off/ratios.txt", quote="\"")




data_hm <- c(hm4$V1)      
data_ht <- c(ht4$V1)
ratio <- c(ratio$V1)
length <- 20000000
 
d1 <-density(data_hm, adjust = 2, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 2, kernel = c("gaussian"))
d3 <- density(ratio, adjust = 2, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "  ", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = 'brown3') ##Homozygous 
lines(d2, col = "royalblue2", lty=2)
lines(d3, col = "gray46", lwd =5)

legend("topright",col=c('brown3', "royalblue2", "grey46"),lwd=1,lty=1:2, cex=1.2, pt.cex = 2,
       legend=c("Homozygous SNP density","Heterozygous SNP density", "Ratio Hom/Het"), bty="n")




