ratio_positions <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/BCF2_3_v2/ratio_positions.csv")
library("ggplot2", lib.loc="/usr/local/Cellar/r/3.1.1/R.framework/Versions/3.1/Resources/library")

x <- ratio_positions$Position
y <- ratio_positions$Ratio

normalized = (y-min(y))/(max(y)-min(y))


g <- ggplot(ratio_positions, aes(x = Position, normalized)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab("Chromosome 4 position (bp)") + ylab("Normalised ratio")
g + theme_bw()
options(scipen=10)

