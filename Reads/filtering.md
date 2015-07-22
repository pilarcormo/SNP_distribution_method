SNP density after parental filtering and centromere removal
====

The homozygous and heterozygous SNP densities obtained after parental filtering and centromere removal were plotted toegether with the ratio signal to identify the high density peaks surrounding the causative mutation position. The R code used to plot the densities is available at [https://github.com/pilarcormo/SNP_distribution_method/blob/master/R_scripts/Overlap_densities.R](https://github.com/pilarcormo/SNP_distribution_method/blob/master/R_scripts/Overlap_densities.R):

```
d1 <-density(data_hm, adjust = 2, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 2, kernel = c("gaussian"))
d3 <- density(ratio, adjust = 2, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "  ", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = 'magenta') ##Homozygous 
lines(d2, col = "royalblue2", lty=2)
lines(d3, col = "gray46", lwd =5)
legend("topright",col=c('magenta', "royalblue2", "grey46"),lwd=1,lty=1:2, cex=1.2, pt.cex = 2,
       legend=c("Homozygous SNP density","Heterozygous SNP density", "Ratio Hom/Het"), bty="n")
 ```


### BCF2 

[Allen et al](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772335/#SM3)

EMS mutagenesis. BCF2 reads and Landsberg erecta (Ler) background. 

**Chromosome 3**

```
hm_bcf2 <- read.table("~/SNP_distribution_method/BCF2/BCF2_chromosome3/interesting_3/hm_nocen.txt", quote="\"")
ht_bcf2 <- read.table("~/NP_distribution_method/BCF2/BCF2_chromosome3/interesting_3/ht_nocen.txt", quote="\"")
ratio_bcf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_bcf2$V1)      
data_ht <- c(ht_bcf2$V1)
ratio <- c(ratio_bcf2$V1)
length <- 23459830 #chromosome 3 length
```
![Image](BCF2/BCF2_chromosome3/Interesting_3/Rplot.withratio.png)

###OCF2 

[Galvão et al](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)

EMS mutagenesis. OCF2 reads and  mir159a parent as background. 

**Chromosome 2**
```
hm_ocf2 <- read.table("~/SNP_distribution_method/BCF2/OCF2_chromosome2/Interesting_2/hm_nocen.txt", quote="\"")
ht_ocf2 <- read.table("~/NP_distribution_method/BCF2/OCF2_chromosome2/Interesting_2/ht_nocen.txt", quote="\"")
ratio_ocf2 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_ocf2$V1)      
data_ht <- c(ht_ocf2$V1)
ratio <- c(ratio_ocf2$V1)
length <- 19698289 #chromosome 2 length
```

![Image](OCF2/OCF2_chromosome2/Interesting_2/Rplot.ratio_nocen.png)

###sup1

[Uchida et al](http://pcp.oxfordjournals.org/content/52/4/716.long)

EMS mutagenesis. sup#1 and sup#2 mutants. Arabidopsis Wassilewskija (Ws) and Col-Tasaka (Col-T) background. 

**Chromosome 4**
```
hm_sup1 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/hm_nocen.txt", quote="\"")
ht_sup1 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/ht_nocen.txt", quote="\"")
ratio_sup1 <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb/SDM_0/ratios.txt", quote="\"")
)
data_hm <- c(hm_sup1$V1)      
data_ht <- c(ht_sup1$V1)
ratio <- c(ratio_sup1$V1)
length <- 18585056 #chromosome 4 length
```

![Image](Aw_sup1-2/Variant_calling/sup1_2_4/Rplot.withratio.png)

###bak1-5 mutants

#####mob1
**Chromosome 5**

```
hm_b <- read.table("~/SNP_distribution_method/m_mutants/B_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
ht_b <- read.table("~/NP_distribution_method/m_mutants/B_chromosome5/interesting_5/ht_nocen.txt", quote="\"")
ratio_b <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_b$V1)      
data_ht <- c(ht_b$V1)
ratio <- c(ratio_b$V1)
length <- 26975502 #chromosome 5 length
```

![Image](m_mutants/B_chromosome5/interesting_5/Rplot.withratio.png)

#####mob2
**Chromosome 5**

```
hm_c <- read.table("~/SNP_distribution_method/m_mutants/C_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
ht_c <- read.table("~/NP_distribution_method/m_mutants/C_chromosome5/interesting_5/ht_nocen.txt", quote="\"")
ratio_c <- read.table("~/SNP_distribution_method/arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb/SDM_0/ratios.txt", quote="\"")
data_hm <- c(hm_c$V1)      
data_ht <- c(ht_c$V1)
ratio <- c(ratio_c$V1)
length <- 26975502 #chromosome 5 length
```

![Image](m_mutants/C_chromosome5/interesting_5/Rplot.withratio.png)






