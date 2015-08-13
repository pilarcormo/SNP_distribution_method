library(ggplot2)
library(grid)
library(gridExtra)
library(ggthemes)
library(RColorBrewer)

density_sum_back <- read.csv("~/SNP_distribution_method/Reads/density_sum_back.csv")
density_sum_out <- read.csv("~/SNP_distribution_method/Reads/density_sum_out.csv")

density_sum_back$cross <- factor(density_sum_back$cross, labels = c("BCF2", "mob1", "mob2"))
density_sum_out$cross <- factor(density_sum_out$cross, labels = c( "OCF2", "sup1"))

#density_sum_out$type <- factor(density_sum_out$type, levels = c("Raw VCF","Parental filtering","Centromere removal", "Candidate SNPs"))
#density_sum_back$type <- factor(density_sum_back$type, levels = c("Raw VCF","Parental filtering","Centromere removal", "Candidate SNPs"))


Palette <- brewer.pal(2,"Set1")
Palette2<-brewer.pal(3,"Set2")

g <- ggplot(density_sum_back, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_back, cross %in% c("BCF2", "mob1", "mob2")), aes(group = cross)) + scale_x_discrete(limits=c("Raw VCF","Parental filtering","Centromere removal"))  +labs(x = "",  y = "Absolute SNP density")  
g <- g +theme_bw(base_size = 20)
g <- g + scale_color_manual(values=Palette2, name  ="Back-cross") + scale_shape_discrete(name  ="Back-cross")+  theme(legend.position="none")


t <- ggplot(density_sum_back, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_back, cross %in% c("BCF2", "mob1", "mob2")), aes(group = cross)) + scale_x_discrete(limits=c("Centromere removal", "Candidate SNPs"))+ ylim(0, 100)  + labs(x = "",  y = "Absolute SNP density") 
t <- t + theme_bw(base_size = 20) 
t <- t +  scale_color_manual(values=Palette2, name  ="Back-cross") + scale_shape_discrete(name  ="Back-cross") 



s <- ggplot(density_sum_out, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_out, cross %in% c("OCF2", "sup1")), aes(group = cross)) + scale_x_discrete(limits=c("Centromere removal", "Candidate SNPs"))+ ylim(0, 10000) + labs(x = " ",  y = "Absolute SNP density")
s <- s +theme_bw(base_size = 20) 
s <- s +  scale_color_manual(values=Palette, name  ="Out-cross")+ scale_shape_discrete(name  ="Out-cross")


h <- ggplot(density_sum_out, aes(x = type, y = density, shape = cross, colour = cross)) + geom_point(size=3) + geom_line(data = subset(density_sum_out, cross %in% c("OCF2", "sup1")), aes(group = cross)) + scale_x_discrete(limits=c("Raw VCF","Parental filtering","Centromere removal"))+ labs(x = " ",  y = "Absolute SNP density")  
h <- h + theme_bw(base_size = 20) 
h <- h+ scale_color_manual(values=Palette, name  ="Out-cross")+ scale_shape_discrete(name  ="Out-cross") + theme(legend.position="none")


ts <- grid.arrange(t, s,  ncol=2, as.table =TRUE)
gh <- grid.arrange(g, t, h, s, ncol=2, as.table =TRUE)


