
hm_bcf2 <- read.table("~/SNP_distribution_method/Reads/BCF2/BCF2_chromosome3/interesting_3/hm_nocen.txt", quote="\"")
ht_bcf2 <- read.table("~/SNP_distribution_method/Reads/BCF2/BCF2_chromosome3/interesting_3/ht_nocen.txt", quote="\"")
ratio_bcf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/10kb_contigs/bcf2_nocen_chr3_10kb/autoratio_1//ratios.txt", quote="\"")
data_hm <- c(hm_bcf2$V1)      
data_ht <- c(ht_bcf2$V1)
ratio <- c(ratio_bcf2$V1)
length <- 23459830 #chromosome 3 length

hm_ocf2 <- read.table("~/SNP_distribution_method/Reads/OCF2/OCF2_chromosome2/Interesting_2/hm_nocen.txt", quote="\"")
ht_ocf2 <- read.table("~/SNP_distribution_method/Reads/OCF2/OCF2_chromosome2/Interesting_2/ht_nocen.txt", quote="\"")
ratio_ocf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/10kb_contigs/ocf2_nocen_chr2_10kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_ocf2$V1)      
data_ht <- c(ht_ocf2$V1)
ratio <- c(ratio_ocf2$V1)
length <- 19698289 #chromosome 2 length


hm5 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/hm_nocen.txt", quote="\"")
ht5 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/ht_nocen.txt", quote="\"")
ratio <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/2-5kb_contigs/sup1_nocen_chr4_5kb/sup1_SDM_0/ratios.txt", quote="\"")

options(scipen = 10)

#ratio <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/ratios_BC.txt", quote="\"")

ratio2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/old/sup1_nocen_chr4_100kb/SDM_off/ratios.txt", quote="\"")
data_hm <- c(hm5$V1)      
data_ht <- c(ht5$V1)
ratio <- c(ratio$V1)
ratio2 <- c(ratio$V1)
length <- 18585056


########PLOT DENSITIES######## 

d1 <-density(data_hm, adjust = 1, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 1, kernel = c("gaussian"))
d3 <- density(ratio, adjust = 2, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "Chromosome 4", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = 'magenta') ##Homozygous 
lines(d2, col = "royalblue2", lty=2)
lines(d3, col = "gray46", lwd =5)

legend("topright",col=c('magenta', "royalblue2", "grey46"),lwd=1,lty=1:2, cex=1.2, pt.cex = 2,
       legend=c("Homozygous SNP density","Heterozygous SNP density", "Ratio Hom/Het"), bty="n")
################ 

