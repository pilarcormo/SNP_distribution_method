library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)

deviations <- read.csv("~/SNP_distribution_method/Small_genomes/1-15Mb.csv")
deviations_30 <- read.csv("~/SNP_distribution_method/Small_genomes/30Mb.csv")

deviations$Genome_size <- factor(deviations$Genome_size, labels = c("1", "3", "5", "7", "9", "11", "13", "15"))
deviations$Contigs <- factor(deviations$Contigs, labels = c("1300", "700"))
deviations_30$Genome_size <- factor(deviations_30$Genome_size, labels = c("30"))
deviations_30$Contigs <- factor(deviations_30$Contigs, labels = c("4000", "2000", "1000"))

Palette <- c('darkblue',"royalblue2")
Palette2 <- c('darkgreen',"green4", "green3")

g <- ggplot(deviations, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) + geom_jitter(size=3) + ylim(0, 3.2) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette) + theme(text = element_text(size=20))
h <- ggplot(deviations_30, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) + geom_jitter(size=3) + ylim(0, 3.2) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette2) + theme(text = element_text(size=20))
gh <- grid.arrange(g , h, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)

#####Theme variations
##Economist
g <- g + theme_economist(base_size = 25) 
h <- h + theme_economist(base_size = 25) 
##Wsj
g <- g + theme_wsj()+ scale_colour_wsj("colors6")
h <- h + theme_wsj()+ scale_colour_wsj("colors6")
##Stata
g <- g + theme_stata(base_size = 20) + scale_color_stata() 
h <- h + theme_stata(base_size = 20) + scale_color_stata() 

gh <- grid.arrange(g , h, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)
