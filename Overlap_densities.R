d1 <- density(hm1$V1)
d2 <- density(ht1$V1)
plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", xlab = "Chromosome length", 
     ylab = "Density")
lines(d1, col = "red")
lines(d2, col = "blue")

