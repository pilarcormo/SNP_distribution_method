library(ggplot2)
options(scipen = 10)
library(RColorBrewer)
library(grid)
library(gridExtra)
######DATA####
contigs <- read.csv("~/SNP_distribution_method/Contigs/contigs.csv")
summary(contigs)
illumina <- read.csv("~/SNP_distribution_method/Contigs/illumina_hiseq.csv")
summary(illumina)
non_illumina <- read.csv("~/SNP_distribution_method/Contigs/non_illumina.csv")
summary(non_illumina)
########

##Effect of technology used on coverage and N50 contig 
palette2 <- brewer.pal(9,"Set1")
pal <- colorRampPalette(palette2)
coverage <- ggplot(contigs, aes(x = Coverage, y = N50_contig, colour = Technology_used)) + geom_point(size = 3) + scale_colour_manual(values=palette2) +labs(x = "Coverage", y = "contig N50") + theme_bw() 
coverage
color_technology <- ggplot(contigs, aes(x = Sequence_lengths, y = N50_contig, colour = Technology_used)) + geom_point(size = 3)  + scale_colour_manual(values=palette2, name = "Technology")  + labs(x = "Genome size (bp)", y = "N50 contig") + theme_bw(base_size = 15)
######

##Probability of N50 contig vs Genome size 
genome_length <- 1000000000
#density <- density(contigs$N50_contig)$y
density <- density(illumina$N50_contig)$y
genome <- (1:512)*(genome_length/512)
df <- data.frame(density, genome)
density_plot <- ggplot(df, aes(x = genome, y = density))  + geom_point(colour ='#F781BF' ) + labs(x = "Genome size (bp)", y = "density of N50") + theme_bw()
density_plot <- density_plot + theme_minimal(base_size = 15)
########         

            
########log(N50) vs log(Genome size) LINEAR REGRESSION and GAM ##########

leg_r2 <- function(k)
{
  legend(x = "topleft", bty = "n",
         legend = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                             list(a = format(coef(k)[1], digits = 2), 
                                  b = format(coef(k)[2], digits = 2), 
                                  r2 = format(summary(k)$r.squared, digits = 3))))
  as.character(as.expression(legend));
}

x_illumina <- log(illumina[illumina$Sequence_lengths > 20000000,]$Sequence_lengths)
y_illumina <- log(illumina[illumina$Sequence_lengths > 20000000,]$N50_contig)
df_logs <- data.frame(x_illumina, y_illumina)

lm_eqn(df_logs)

lm(illumina_x ~ illumina_y)
summary(lm(illumina_x ~ illumina_y, data = df_logs))$r.squared


x_bla <- df_logs[df_logs$x_illumina < 20.5,]$x_illumina
y_bla <- df_logs[df_logs$x_illumina < 20.5,]$y_illumina
df_logs2 <- data.frame(x_bla, y_bla)
fit <- lm(x_bla ~ y_bla, data = df_logs2)
len <- leg_r2(fit)

p <- function(df_logs2){
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(lm(y_bla ~ x_bla))[1], digits = 2), b = format(coef(lm(y_bla ~ x_bla))[2], digits = 2), r2 = format(summary(lm(y_bla ~ x_bla, data = df_logs2))$r.squared, digits = 3)))
as.character(as.expression(eq))
}

p2 <- function(df_logs){
  eq <- substitute(italic(r)^2~"="~r2, list(r2 = 0.807))
  as.character(as.expression(eq))
}

s <- ggplot(df_logs, aes(x = x_illumina, y = y_illumina))  + geom_jitter(size = 3, colour = "#F781BF") +labs(x = "ln Genome size (bp)", y = "ln N50 contig") + theme_bw() 
s <- s + geom_smooth(method = "lm", data=subset(df_logs, x_illumina < 20.5), se= FALSE, size = 1, colour = "black") + annotate("text", x = 19, y = 11, label = p(df_logs2), parse = TRUE, colour = "black")

s + geom_smooth(method = "gam", formula = y ~ s(x), size = 1, colour = "#F781BF", fill = "gainsboro" ) + annotate("text", x = 21.8, y = 11.5, label = p2(df_logs)  , parse = TRUE, colour = "#F781BF")

########


############GAM model##############
x <- log(illumina[illumina$Sequence_lengths > 20000000,]$Sequence_lengths)

y <- log(illumina[illumina$Sequence_lengths > 20000000,]$N50_contig)
df <- data.frame(x, y)
library(mgcv)
model <- gam(y ~ s(x))
summary(model)
coef(gam(y ~ s(x)))
fit <- lm(x_bla ~ y_bla, data = df_logs2)

##########        

#########N50 Densities, mean and median########
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


