QQ-plots 
====
R code plot the probability plots (QQ-plots) with the linear regression can be found at [https://github.com/pilarcormo/SNP_distribution_method/blob/master/R_scripts/qqplot.R](https://github.com/pilarcormo/SNP_distribution_method/blob/master/R_scripts/qqplot.R)

```
##qqplot and qqline2
qqline2 <- function(x, y, probs = c(0.25, 0.75), qtype = 7, ...)
{
  stopifnot(length(probs) == 2)
  x2 <- quantile(x, probs, names=FALSE, type=qtype, na.rm = TRUE)
  y2 <- quantile(y, probs, names=FALSE, type=qtype, na.rm = TRUE)
  slope <- diff(y2)/diff(x2)
  int <- y2[1L] - slope*x2[1L]
  abline(int, slope, ...)
}
leg_r2 <- function(k) #Generate legend with the equation of the line and the value of r square
{
  legend(x = "topleft", bty = "n",
         legend = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                             list(a = format(coef(k)[1], digits = 2), 
                                  b = format(coef(k)[2], digits = 2), 
                                  r2 = format(summary(k)$r.squared, digits = 3))))
}
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
```

###1. sup1 
[Uchida et al](http://pcp.oxfordjournals.org/content/52/4/716.long)


#####Chromosome 4
```
hm_sup1 <- read.table("~/SNP_distribution_method/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/hm_nocen.txt", quote="\"")
y5 <- c(hm_sup1$V1)
q5 <- qqplot_line(y5, "sup#1")
```

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/Rplot.qqsup1.png?raw=true)


###OCF2
[Galvão et al](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)
 
 
##### Chromosome 2

```
hm_OCF2 <- read.table("~/SNP_distribution_method/Reads/OCF2/chromosome2/Interesting_2/hm_nocen.txt", quote="\"")
y4 <- c(hm_OCF2$V1)
q4 <- qqplot_line(y4, "OCF2")
```

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/OCF2/OCF2_chromosome2/Interesting_2/Rplot.qqof.png?raw=true)


###mob mutants 


#####B - Chromosome 5
```
hm_B <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
y1 <- c(hm_B$V1)
q1 <- qqplot_line(y1, "mob 1")
```
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/B_chromosome5/interesting_5/Rplot.qqm1.png?raw=true)


#####C - Chromosome 5
```
hm_C <- read.table("~/SNP_distribution_method/Reads/m_mutants/C_chromosome5/interesting_5/hm_nocen.txt", quote="\"")
y2 <- c(hm_C$V1)
q2 <- qqplot_line(y2, "mob 2")
```

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/C_chromosome5/interesting_5/Rplot.qqm2.png?raw=true)

###BCF2

##### Chromosome 3
```
hm_BCF2 <- read.table("~/SNP_distribution_method/Reads/BCF2/chromosome3/Interesting_3/hm_nocen.txt", quote="\"")
y3 <- c(hm_BCF2$V1)
q3 <- qqplot_line(y3, "BCF2")
```

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/BCF2/BCF2_chromosome3/Interesting_3/Rplot.qqbc.png?raw=true)


