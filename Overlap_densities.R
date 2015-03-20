data_hm <- c(hm5$V1)
data_ht <- c(ht5$V1)
length <- 27000000

d1 <-density(data_hm, adjust = 3, kernel = c("gaussian"))
d2 <- density(data_ht, adjust = 5, kernel = c("gaussian"))
p1 <- plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", main = "Chromosome 5")
lines(d1, col = "magenta2") ##Homozygous 
lines(d2, col = "royalblue2", lty=2)

legend("topright",col=c("magenta2", "royalblue2"),lwd=1,lty=1:3,legend=c("Homozygous SNP density","Heterozygous SNP density"), bty="n")




