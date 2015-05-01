#encoding: utf-8

require_relative 'write_it'
require 'rinruby'

# hm = "Arabidopsis_datasets/No_centromere/sup1_nocen_chr4_100kb/global_perm_hm0.txt"
# hm2 = "Reads/Aw_sup1-2/Variant_calling/sup1_2_4/hm_nocen.txt"


class Plot
  def self.exp_vs_hyp_densities(expected, hypothetical, location)
    myr = RinRuby.new(echo = false)
    myr.assign "hypothetical", hypothetical
    myr.assign "expected", expected
    myr.assign "location", location
    myr.eval 'png(paste(location,"/", "qqplot_exp_hyp",".png", sep=""))
    qqline2 <- function(x, y, probs = c(0.25, 0.75), qtype = 7, ...)
      {
        stopifnot(length(probs) == 2)
        x2 <- quantile(x, probs, names=FALSE, type=qtype, na.rm = TRUE)
        y2 <- quantile(y, probs, names=FALSE, type=qtype, na.rm = TRUE)
        slope <- diff(y2)/diff(x2)
        int <- y2[1L] - slope*x2[1L]
        abline(int, slope, ...)
      }

    leg_r2 <- function(k)
      {
        legend(x = "topleft", bty = "n",
               legend = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                   list(a = format(coef(k)[1], digits = 2), 
                                        b = format(coef(k)[2], digits = 2), 
                                        r2 = format(summary(k)$r.squared, digits = 3))))
      }

    y <- expected
    x <- hypothetical
    df <- data.frame(x, y)
    V = qqplot(x, y, main="Q-Q Plot", ylab="Expected homozygous SNP density", xlab = "hypothetical homozygous SNP density after SDM")
    l <- qqline2(x, y, col = 6) 
    fg <- data.frame(V$x, V$y)
    k <- lm(V$y ~ V$x)
    len <- leg_r2(k)'
    myr.quit
  end 

  def self.densities(hm, ht, ratio, length, location)
    myr = RinRuby.new(echo = false)
    myr.hm =  hm
    myr.ht =  ht
    myr.ratio = ratio
    myr.length = length
    myr.eval 'png(paste(location,"/", "qqplot_exp_hyp",".png", sep="")
    d1 <-density(hm, adjust = 1, kernel = c("gaussian"))
    d2 <- density(ht, adjust = 1, kernel = c("gaussian"))
    d3 <- density(ratio, adjust =1 , kernel = c("gaussian"))
    p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "Densities", xlim =c(0,length), xlab = " ", ylab = " ")
    lines(d1, col = "magenta2") ##Homozygous 
    lines(d2, col = "royalblue2", lty=2)
    lines(d3, col = "gray46", lwd =5)
    legend("topright",col=c("magenta2", "royalblue2", "grey46"),lwd=1,lty=1:2,
       legend=c("Expected homozygous SNP density","Expected heterozygous SNP density", "Expected hom/het ratio"), bty="n")'
    myr.quit
  end
end 

# reads = "Aw_sup1-2/Variant_calling"
# hypothetical = WriteIt.file_to_ints_array("../arabidopsis_datasets/No_centromere/sup1_nocen_chr4_100kb/global_perm_hm0.txt") # Get SNP distributions
# expected = WriteIt.file_to_ints_array("../Reads/Aw_sup1-2/Variant_calling/sup1_2_4/hm_nocen.txt")
# location = "../Reads/#{reads}"
# Plot::exp_vs_hyp_densities(expected, hypothetical, location)
