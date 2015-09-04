library(ggplot2)
options(scipen = 10)
library(RColorBrewer)
library(grid)
library(gridExtra)

contigs <- read.csv("~/SNP_distribution_method/Contigs/contigs.csv")
summary(contigs)
illumina <- read.csv("~/SNP_distribution_method/Contigs/illumina_hiseq.csv")
summary(illumina)
non_illumina <- read.csv("~/SNP_distribution_method/Contigs/non_illumina.csv")
summary(non_illumina)

##Effect of technology used on coverage and N50 contig 
palette2 <- brewer.pal(9,"Set1")
pal <- colorRampPalette(palette2)
coverage <- ggplot(contigs, aes(x = Coverage, y = N50_contig, colour = Technology_used)) + geom_point(size = 3) + scale_colour_manual(values=palette2) +labs(x = "Coverage", y = "contig N50") + theme_bw() 
coverage
color_technology <- ggplot(contigs, aes(x = Sequence_lengths, y = N50_contig, colour = Technology_used)) + geom_point(size = 3)  + scale_colour_manual(values=palette2, name = "Technology")  + labs(x = "Genome size (bp)", y = "N50 contig") + theme_bw(base_size = 15)

##Probability of N50 contig vs Genome size 
genome_length <- 1000000000
#density <- density(contigs$N50_contig)$y
density <- density(illumina$N50_contig)$y
genome <- (1:512)*(genome_length/512)
df <- data.frame(density, genome)
density_plot <- ggplot(df, aes(x = genome, y = density))  + geom_point(colour ='#F781BF' ) + labs(x = "Genome size (bp)", y = "density of N50") + theme_bw()
density_plot <- density_plot + theme_minimal(base_size = 15)
####           

            
##log(N50) vs log(Genome size) 
x <- log(illumina[illumina$Sequence_lengths > 20000000,]$Sequence_lengths)
y <- log(illumina[illumina$Sequence_lengths > 20000000,]$N50_contig)
colour <- contigs[contigs$Sequence_lengths > 20000000,]$Technology_used

df_logs <- data.frame(x, y)
r <- ggplot(df_logs, aes(x = x, y = y))  + geom_jitter(size = 3, colour = "#F781BF") +labs(x = "log Genome size (bp)", y = "log N50 contig") + theme_bw() 

r <- r + geom_smooth(colour = "#F781BF", fill = "grey85") + theme_minimal(base_size = 15)

##GAM model
x <- log(illumina[illumina$Sequence_lengths > 20000000,]$Sequence_lengths)

y <- log(illumina[illumina$Sequence_lengths > 20000000,]$N50_contig)
df <- data.frame(x, y)
library(mgcv)
model <- gam(y ~ s(x))
summary(model)
coef(gam(y ~ s(x)))

plot(model, residuals=T, pch=19,
     scheme=1, col='#F781BF', shade=T, xlab ="log (Genome size (bp))", ylab = "log (N50 contig)")
legend(x = "topright", lwd=1,
       c("R-sq = 0.807"),
       col = c("#F781BF"), bty="n")
                 
##N50 Densities, mean and median
den_illu <- density(illumina$N50_contig)
den_con <- density(contigs$N50_contig)
p1 <- plot(range(den_con$x, den_illu$x), range(den_con$y, den_illu$y), type = "n", main = "N50 length distribution", xlab = "N50 contig length", ylab = " ")
lines(den_con, col = "#377EB8")
lines(den_illu, col = "#F781BF") 

abline(v = median(illumina$N50_contig),
       col = "#F781BF",
       lty=2)
abline(v = median(contigs$N50_contig), 
       col = "#377EB8", 
       lty=2)

legend(x = "topright", lwd=1,
       c("Non-Illumina", "Illumina HiSeq"),
       col = c("#377EB8", "#F781BF"),
       lty = c(2, 2), bty="n")

legend("topright",col=c("#377EB8", "#F781BF"),lwd=1,
       legend=c("Other", "Illumina HiSeq"), bty="n")
#########


