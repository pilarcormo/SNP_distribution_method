library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)

deviations <- read.csv("~/SNP_distribution_method/Small_genomes/arabidopsis_datasets/1-15Mb.csv")
deviations_30 <- read.csv("~/SNP_distribution_method/Small_genomes/arabidopsis_datasets/30Mb.csv")

deviations$Genome_size <- factor(deviations$Genome_size, labels = c("1", "3", "5", "7", "9", "11", "13", "15"))
deviations$Contigs <- factor(deviations$Contigs, labels = c("1300", "700"))
deviations_30$Genome_size <- factor(deviations_30$Genome_size, labels = c("30"))
deviations_30$Contigs <- factor(deviations_30$Contigs, labels = c("4000", "2000", "1000"))

Palette <- c('red',"royalblue2")
Palette2 <- c('darkgreen',"green4", "green3")

p15 <- ggplot(deviations, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) + geom_jitter(size=3) + ylim(0, 3.2) + labs(x = "Genome length (Mb)",  y = "Shift from 'real' SNP as % of genome length") + theme_bw() + scale_colour_manual(values=Palette, name  ="Number of contigs") + scale_shape_discrete(name  ="Number of contigs") + theme(text = element_text(size=20))
p30 <- ggplot(deviations_30, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) + geom_jitter(size=3) + ylim(0, 3.2) + labs(x = "Genome length (Mb)",  y = "Shift from 'real' SNP as % of genome length") + theme_bw() + scale_colour_manual(values=Palette2, name ="Number of contigs") + scale_shape_discrete(name  ="Number of contigs") + theme(text = element_text(size=20))


grid <- grid.arrange(p15 , p30, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)

#####Theme variations
##Economist
g <- g + theme_economist(base_size = 25) 
h <- h + theme_economist(base_size = 25) 
##Wsj
p15 <- p15 + theme_minimal(base_size = 15)
p30 <- p30 + theme_minimal(base_size = 15)
##Stata
g <- g + theme_stata(base_size = 15) + scale_color_stata() 
h <- h + theme_stata(base_size = 15) 

grid <- grid.arrange(p15 , p30, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)
