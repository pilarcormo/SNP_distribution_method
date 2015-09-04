##Option 1- qqplot and qqline2
qqline2 <- function(x, y, probs = c(0.25, 0.75), qtype = 7, ...)
{
  stopifnot(length(probs) == 2)
  x2 <- quantile(x, probs, names=FALSE, type=qtype, na.rm = TRUE)
  y2 <- quantile(y, probs, names=FALSE, type=qtype, na.rm = TRUE)
  slope <- diff(y2)/diff(x2)
  int <- y2[1L] - slope*x2[1L]
  abline(int, slope, ...)
}
##Generate legend with the equation of the line and the value of r square
leg_r2 <- function(k)
{
  legend(x = "topleft", bty = "n",
         legend = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                             list(a = format(coef(k)[1], digits = 2), 
                                  b = format(coef(k)[2], digits = 2), 
                                  r2 = format(summary(k)$r.squared, digits = 3))))
}

hm <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
hm2 <- read.table("~/SNP_distribution_method/Reads/m_mutants/B_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
hm3 <- read.table("~/SNP_distribution_method/Reads/BCF2/Interesting_3/hm_nocen.txt", quote="\"")

hm_pa <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/sup1_chromosome4/hm.txt", quote="\"")

hm_cen <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/hm.txt", quote="\"")

hm_pre <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/filter2_chromosome4/hm.txt", quote="\"")

options(scipen = 10)
y1 <- c(hm$V1)
y2 <- c(hm2$V1)
y3 <- c(hm3$V1)
y4 <- c(hm4$V1)
y5 <- c(hm5$V1)

y_pa <- c(hm_pa$V1)
y_cen <- c(hm_cen$V1)
y_pre <- c(hm_pre$V1)

qqplot_line <- function(y, title)
{
  x <- rnorm(length(y), mean(y), sd(y))
  df <- data.frame(x, y)
  V = qqplot(x, y, main = title, ylab="Real homozygous SNP density", xlab = "Theoretical normal disitribution")
  l <- qqline2(x, y, col = 'brown3') 
  fg <- data.frame(V$x, V$y)
  k <- lm(V$y ~ V$x)
  len <- leg_r2(k)
} 

q1 <- list(qqplot_line(y_pre, "sup1 pre-filtering"))
q2 <- list(qqplot_line(y_pa, "sup1 after parental filtering"))
q3 <- list(qqplot_line(y_cen, "sup1 after centromere removal"))

q3 <- qqplot_line(y3, "BCF2")
q4 <- qqplot_line(y4, "OCF2")

q5 <- qqplot_line(y5, "sup#1")


grid.arrange(arrangeGrob(q1, q2), ncol = 2)

##Standard deviation, kurtosis and skewness of the distribution
library(moments)
sd(y)
kurtosis(y)
skewness(y)

##QQnorm
n <- qqnorm(y); qqline(y, col = 6)   

##Option 2- ggplot2, stat_qq and geom_line

ggQQ <- function(df) # argument: a linear model
{
  y <- quantile(df$y, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  p <- ggplot(df, aes(sample=y)) +
    stat_qq(alpha = 0.5, distribution = nbinom) +
    geom_abline(slope = slope, intercept = int, color="blue") + theme_bw()
}

p <- ggQQ(df)

