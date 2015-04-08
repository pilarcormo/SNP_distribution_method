#ratio_positions <- read.csv("~/SNP_distribution_method/arabidopsis_datasets/Genomes_SDM/chr1_A_1/ratio_positions0.csv")
x <- ratio_positions0$Position
y <- ratio_positions0$Ratio

normalized = (y-min(y))/(max(y)-min(y))

g <- ggplot(ratio_positions0, aes(Position, Ratio)) + geom_point(colour="darkblue", shape=21, size = 4) + xlab("Genome position (bp)") + ylab("Normalised ratio") + xlim(0, 30000000)
g + theme_bw()

