#encoding: utf-8

require_relative 'write_it'
require 'rinruby'

class Plot
  def self.qqplot(hypothetical, location, title, ylabel, xlabel)
    myr = RinRuby.new(echo = false)
    myr.assign "hypothetical", hypothetical
    myr.assign 'title', title
    myr.assign 'xlabel', xlabel
    myr.assign 'ylabel', ylabel
    myr.assign "location", location
    myr.eval 'png(paste(location,"/", "qqplot_exp_hyp",".png", sep=""), width=500,height=500)
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

    x <- hypothetical
    y <- rnorm(length(x), mean(x), 5000000)
    df <- data.frame(x, y)
    V = qqplot(x, y, main="Q-Q Plot", ylab=ylabel, xlab =xlabel)
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
    myr.location = location 
    myr.eval 'png(paste(location,"/", "Hypothetical densities",".png", sep=""), width=800,height=500)
    options(scipen = 10) 
    d1 <-density(hm, adjust = 1, kernel = c("gaussian"))
    d2 <- density(ht, adjust = 1, kernel = c("gaussian"))
    d3 <- density(ratio, adjust =0.25 , kernel = c("gaussian"))
    p1 <- plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n", main = "Densities", xlim =c(0,length), xlab = " ", ylab = " ")
    lines(d1, col = "magenta2") ##Homozygous 
    lines(d2, col = "royalblue2", lty=2)
    lines(d3, col = "gray46", lwd =5)
    legend("topright",col=c("magenta2", "royalblue2", "grey46"),lwd=1,lty=1:2,
        legend=c("Homozygous SNP density","Heterozygous SNP density", "hom/het ratio"), bty="n")'
    myr.quit
  end

  def self.get_ylim(density, genome_length)
      myr = RinRuby.new(echo = false)
      myr.assign 'density', density
      myr.assign 'genome_length', genome_length
      myr.eval 'plot((1:512)*(genome_length/512), density(density)$y)'
      ylim = myr.pull 'par("yaxp")[2] + par("yaxp")[2]/par("yaxp")[3]'
      myr.quit
  return ylim
  end 

  def self.comparison(expected, hypothetical, length, location, ylim)
    myr = RinRuby.new(echo = false)
    myr.hypothetical = hypothetical
    myr.expected =  expected
    myr.length = length
    myr.location = location 
    myr.ylim = ylim
    myr.eval 'png(paste(location,"/", "ratios",".png", sep=""), width=500,height=500)
    options(scipen = 10) 
    d1 <-density(hypothetical, adjust = 1, kernel = c("gaussian"))
    d2 <-density(expected, adjust = 1, kernel = c("gaussian"))
    p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = "Ratios", xlim =c(0,length), xlab = " ", ylab = " ")
    lines(d1, col = "slategray4", lwd =3)  
    lines(d2, col = "steelblue3", lwd =3, lty=2)
    legend("topright",col=c("slategray4", "steelblue3"),lwd=3, lty=1:2,
        legend=c("Hypothetical ratio","Expected ratio"), bty="n")'
    myr.quit
  end


end 

