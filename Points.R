ratio_positions0 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/OC_chr2_1kb/ratio_positions0_1.csv")
ratio_positions01 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/OC_chr2_1kb/ratio_positions0_0.1.csv")
ratio_positions001 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/OC_chr2_1kb/ratio_positions0_0.01.csv")

ratio_positions2 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/OC_chr2_10kb/ratio_positions0_0.01.csv")

#x <- ratio_positions01$Position
#y <- ratio_positions01$Ratio

#normalized = (y-min(y))/(max(y)-min(y))

ratio1 <- ggplot(ratio_positions2, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab(" ") + ylab("Ratio") + xlim(0, 20000000) + theme_bw()
#ratio01 <- ggplot(ratio_positions01, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab(" ") + ylab("Ratio") + xlim(0, 25000000) + theme_bw() + ggtitle("Factor = 0.1")
#ratio001 <- ggplot(ratio_positions001, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab("Chromosome position (bp)") + ylab("Ratio") + xlim(0, 25000000) + theme_bw() + ggtitle("Factor = 0.01")

#together <- grid.arrange(ratio1, ratio01, ratio001, ncol = 1, main = "1kb contigs")

