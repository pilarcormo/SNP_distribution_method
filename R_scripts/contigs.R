library(ggplot2)
options(scipen = 10)


contigs <- read.csv("~/SNP_distribution_method/Contigs/contigs.csv")

summary(contigs)

illumina <- read.csv("~/SNP_distribution_method/Contigs/illumina_hiseq.csv")

summary(illumina)

non_illumina <- read.csv("~/SNP_distribution_method/Contigs/non_illumina.csv")

summary(non_illumina)

Palette <- c('magenta1','royalblue2',"green3", "green4", "yellowgreen", 'black', 'darkorchid3', 'deepskyblue1', 'gold', 'firebrick1', 'grey', "orange", "skyblue4", "tomato4", 'yellow', 'aquamarine')

g <- ggplot(contigs, aes(x = Coverage, y = N50_contig, colour = Technology_used)) + geom_point(size = 3) + xlim(0, 100) + scale_colour_manual(values=Palette) +labs(x = "Coverage", y = "contig N5") + theme_bw() 

h <- ggplot(contigs, aes(x = Sequence_lengths, y = N50_contig, colour = Technology_used)) + geom_jitter(size = 3) + scale_colour_manual(values=Palette) +labs(x = "Genome size (bp)", y = "N50 contig") + theme_bw()

gh <- grid.arrange(arrangeGrob(g + theme(legend.position="none"), h + theme(legend.position="center")))

                   
den_con <- density(illumina$N50_contig)
den_illu <- density(non_illumina$N50_contig)
p1 <- plot(range(den_con$x), range(den_con$y), type = "n", main = "N50 length distribution", xlab = "N50 contig length", ylab = " ")
lines(den_con, col = "royalblue2")
lines(den_illu, col = "firebrick1") 

abline(v = mean(illumina$N50_contig),
       col = "black",
       lty=2)
abline(v = median(illumina$N50_contig), 
       col = "red", 
       lty=2)

legend(x = "topright", lwd=1,
       c("Mean", "Median"),
       col = c("black", "red"),
       lty = c(2, 2), bty="n")

legend("topright",col=c("firebrick1", "royalblue2"),lwd=1,
       legend=c("Other", "Illumina HiSeq"), bty="n")


median(contigs$N50_contig)
var(contigs$N50_contig)
mean(contigs$N50_contig)

qqnorm(contigs)
qqline(contigs)




