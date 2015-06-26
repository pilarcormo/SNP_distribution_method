library(ggplot2)
library(grid)
library(gridExtra)


deviations <- read.csv("~/Desktop/mutation_1-15Mb.csv")
deviations_30 <- read.csv("~/Desktop/30mb.csv")

deviations$Genome_size <- factor(deviations$Genome_size, labels = c("1", "3", "5", "7", "9", "11", "13", "15"))
deviations$Contigs <- factor(deviations$Contigs, labels = c("1300", "700"))
deviations_30$Genome_size <- factor(deviations_30$Genome_size, labels = c("30"))
deviations_30$Contigs <- factor(deviations_30$Contigs, labels = c("4000", "2000", "1000"))

Palette <- c('brown3',"royalblue2")
Palette2 <- c('green4',"green3", "greenyellow")


g <- ggplot(deviations, aes(x = Genome_size, y = X., shape = Contigs, colour = Contigs)) + geom_point(width = .2) + ylim(0, 1.4) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette) + theme_bw()
h <- ggplot(deviations_30, aes(x = Genome_size, y = X., shape = Contigs, colour = Contigs)) + geom_point(width = .2) + ylim(0, 1.4) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette2) + theme_bw()

gh <- grid.arrange(g , h, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)


