

contigs <- read.csv("~/SNP_distribution_method/contigs.csv")

illumina <- read.csv("~/SNP_distribution_method/illumina_hiseq.csv")


contigs$Technology_used <- factor(Contigs$Technology_used)
contigs$Average_N50 <- factor(deviations$Contigs, labels = c("1300", "700"))

Palette <- c('magenta2','royalblue2',"green3", "green4", "yellowgreen", 'black', 'darkorchid2', 'deepskyblue1', 'gold', 'firebrick1', 'grey', "orange", "skyblue4", "tomato4", 'yellow', 'aquamarine')

g <- ggplot(contigs, aes(x = Sequence_lengths, y = Average_contig_size, colour = Technology_used)) + geom_jitter(size = 3) + scale_colour_manual(values=Palette) +labs(x = "Genome size (Mb)", y = "Average contig size") + theme_bw() 

h <- ggplot(contigs, aes(x = Sequence_lengths, y = N50_contig, colour = Technology_used)) + geom_jitter(size = 3) + scale_colour_manual(values=Palette) +labs(x = "Genome size (Mb)", y = "N50 contig") + theme_bw() + theme(legend.position = "center")

gh <- grid.arrange(arrangeGrob(g + theme(legend.position="none"), h + theme(legend.position="topleft")))


                   
den_con <- density(contigs$Sequence_lengths)
den_illu <- density(illumina$Sequence_lengths)
p1 <- plot(range(den_con$x, den_illu$x), range(den_con$y, den_illu$y), type = "n", main = "N50 length distribution", xlab = "Genome length / N50 contig", ylab = " ")
lines(den_con, col = "royalblue2") 
lines(den_illu, col = "firebrick1") 


legend("topright",col=c("firebrick1", "royalblue2"),lwd=1,
       legend=c("Illumina HiSeq", "All"), bty="n")

p1 <- plot(range(den_con$x), range(den_con$y), type = "n", main = "N50 length distribution", xlab = "N50", ylab = " ")
lines(den_con, col = "firebrick1") 


contigs <- c(contigs$)
qqnorm(contigs)
qqline(contigs)
mean(contigs)
sd(conti)

length(contigs$Average_N50)
mean(contigs$Average_N50)
sd(contigs$Average_N50)

options(scipen = 10)
