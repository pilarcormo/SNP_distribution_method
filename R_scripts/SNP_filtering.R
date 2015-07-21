library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)

density_sum_back <- read.csv("~/SNP_distribution_method/Reads/density_sum_back.csv")
density_sum_out <- read.csv("~/SNP_distribution_method/Reads/density_sum_out.csv")

density_sum_back$cross <- factor(density_sum_back$cross, labels = c("BCF2", "mob1", "mob2"))
density_sum_out$cross <- factor(density_sum_out$cross, labels = c( "OCF2", "sup1"))


g <- ggplot(density_sum_back, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_back, cross %in% c("BCF2", "mob1", "mob2")), aes(group = cross)) + scale_x_discrete(limits=c("Pre-filtering","Parental","Centromere"))+ labs(x = "",  y = "Absolute SNP density") + theme_bw() 
g <- g + scale_colour_discrete(name  ="Back-cross") + scale_shape_discrete(name  ="Back-cross") 
g <- g + theme_stata(base_size = 15)  
h <- ggplot(density_sum_out, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_out, cross %in% c("OCF2", "sup1")), aes(group = cross)) + scale_x_discrete(limits=c("Pre-filtering","Parental","Centromere"))+ labs(x = " ",  y = "Absolute SNP density") + theme_bw() 
h <- h + scale_colour_discrete(name  ="Out-cross") + scale_shape_discrete(name  ="Out-cross") 
h <- h+ theme_stata(base_size = 15) 

gh <- grid.arrange(g , h, ncol=2, as.table =TRUE)

