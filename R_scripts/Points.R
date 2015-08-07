
ratio_positions0_0.1 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/sup1_chr4_10kb/ratio_positions0_0.1.csv")

ratio_positions0_0.01 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/sup1_chr4_10kb/ratio_positions0_0.01.csv")

ratio_positions0_1 <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/Centromere/sup1_chr4_10kb/ratio_positions0_1.csv")

x <- ratio_positions$Position
y <- ratio_positions$Ratio

ratio1 <- ggplot(ratio_positions0_1, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab(" ") + ylab("Ratio") + xlim(5000000, 15000000) + theme_bw() + ggtitle("Factor = 1")
ratio01 <- ggplot(ratio_positions0_0.1, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab(" ") + ylab("Ratio") + xlim(5000000, 15000000) + theme_bw() + ggtitle("Factor = 0.1")
ratio001 <- ggplot(ratio_positions0_0.01, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab(" ") + ylab("Ratio") + xlim(5000000, 15000000) + theme_bw() + ggtitle("Factor = 0.01")



together <- grid.arrange(ratio1, ratio01, ratio001, ratiolen,  ncol = 2, main = "10kb contigs")

