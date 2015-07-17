library(ggthemes)


hyp <- read.table("~/SNP_distribution_method/Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left/Ratio_0_1/hyp_ratios.txt", quote="\"")
exp <- read.table("~/SNP_distribution_method/Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left/Ratio_0_1/ratios.txt", quote="\"")


cex =2 
options(scipen = 10) 
d1 <-density(hyp$V1, adjust = 1, kernel = c("gaussian"))
d2 <-density(exp$V1, adjust = 1, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = " ", xlim =c(0,length), xlab = " ", ylab = " ")
lines(d1, col = "slategray4", lwd =3) 
lines(d2, col = "steelblue3", lwd =3, lty=2)
#legend("topright",col=c("slategray4", "steelblue3"),lwd=3, lty=1:2,
       #legend=c("SDM ratio","Expected ratio"), bty="n", cex=1.3, pt.cex = 2)

par(bg = 'white')
