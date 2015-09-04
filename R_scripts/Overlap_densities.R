
hm_bcf2 <- read.table("~/SNP_distribution_method/Reads/BCF2/BCF2_chromosome3/interesting_3/hm_nocen.txt", quote="\"")
ht_bcf2 <- read.table("~/SNP_distribution_method/Reads/BCF2/BCF2_chromosome3/interesting_3/ht_nocen.txt", quote="\"")
  ratio_bcf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/100kb_contigs/bcf2_nocen_chr3_100kb/Ratio_0_1/ratios.txt", quote="\"")
data_hm <- c(hm_bcf2$V1)      
data_ht <- c(ht_bcf2$V1)
ratio <- c(ratio_bcf2$V1)
length <- 23459830 #chromosome 3 length

hm_ocf2 <- read.table("~/SNP_distribution_method/Reads/OCF2/OCF2_chromosome2/Interesting_2/hm_nocen.txt", quote="\"")
ht_ocf2 <- read.table("~/SNP_distribution_method/Reads/OCF2/OCF2_chromosome2/Interesting_2/ht_nocen.txt", quote="\"")
ratio_ocf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/100kb_contigs/ocf2_nocen_chr2_100kb/SDM_off/ratios.txt", quote="\"")
data_hm <- c(hm_ocf2$V1)      
data_ht <- c(ht_ocf2$V1)
ratio <- c(ratio_ocf2$V1)
length <- 19698289 #chromosome 2 length


hm_sup1 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/hm_nocen.txt", quote="\"")
ht_sup1 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/ht_nocen.txt", quote="\"")
ratio_sup1 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/10kb_contigs/sup1_nocen_chr4_100kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_sup1$V1)      
data_ht <- c(ht_sup1$V1)
ratio <- c(ratio_sup1$V1)
length <- 18585056 #chromosome 4 length

hm_b <- read.table("~/SNP_distribution_method/Reads/m_mutants/B_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
ht_b <- read.table("~/SNP_distribution_method/Reads/m_mutants/B_chromosome5/interesting_5/ht_nocen.txt", quote="\"")
ratio_b <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/100kb_contigs/B_nocen_chr2_100kb/SDM_off/ratios.txt", quote="\"")
data_hm <- c(hm_b$V1)      
data_ht <- c(ht_b$V1)
ratio <- c(ratio_b$V1)
length <- 26975502 #chromosome 5 length

hm_c <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
ht_c <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/ht_nocen.txt", quote="\"")
ratio_c <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/100kb_contigs/C_nocen_chr5_100kb/New_plots_off/ratios.txt", quote="\"")
data_hm <- c(hm_c$V1)      
data_ht <- c(ht_c$V1)
ratio <- c(ratio_c$V1)
length <- 26975502 #chromosome 5 length


########PLOT DENSITIES######## 
options(scipen = 10)
d1 <-density(data_hm, adjust = 1, kernel = c("gaussian"))
d2 <- density(data_ht, adjust =2, kernel = c("gaussian"))
d3 <- density(ratio, adjust = 1, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "mob1", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = 'magenta') ##Homozygous
lines(d2, col = "royalblue2", lty=2)
lines(d3, col = "gray46", lwd =5)

legend("topright",col=c('magenta', "royalblue2", "grey46"),lwd=1,lty=1:2, cex=1.2, pt.cex = 2,
       legend=c("Homozygous SNP density","Heterozygous SNP density", "Ratio Hom/Het"), bty="n")

abline(v = 26450000,
       col = "magenta",
       lty=2)
################ 

