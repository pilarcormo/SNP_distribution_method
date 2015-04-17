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


hm <- read.table("~/SNP_distribution_method/Aw_sup1-2/Variant_calling/sup1_2_1/hm.txt", quote="\"")

x <- rnorm(length(y), mean(y), sd(y))
y <- c(hm$V1)
df <- data.frame(x, y)

V = qqplot(x, y, main="Q-Q Plot for sup#1", ylab="hm SNP density in chromosome 1", xlab = "Theoretical normal disitribution")
l <- qqline2(x, y, col = 6) 
fg <- data.frame(V$x, V$y)
k <- lm(V$y ~ V$x)
len <- leg_r2(k)

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

