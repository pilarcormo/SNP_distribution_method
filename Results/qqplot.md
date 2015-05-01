qqplots I - Understanding the distribution
====

###1. [Uchida et al](http://pcp.oxfordjournals.org/content/52/4/716.long)
#####Chromosome 1 

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_1/Rplot.qq.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_1/Rplot.qnorm.png?raw=true)

#####Chromosome 4

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/Rplot.r2qq.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/Rplot.normal.png?raw=true)


###2. [Galvão et al](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)

Focusing only in the fragments that contribute to the normal distribution

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/OCF2/Interesting_2/Rplot.qqplot.png?raw=true)

###3. m mutants 

#####A - Chromosome 2
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/A_chromosome2/Rplot.qqplot.nocen.png?raw=true)

#####B - Chromosome 5
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/B_chromosome5/Rplot.nocen.qqplot.png?raw=true)

#####C - Chromosome 5
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/C_chromosome5/Rplot.no_cen.qqplot.png?raw=true)
#####F - Chromosome 1
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/F_chromosome1/Rplot.no_cen.qqplot.png?raw=true)


qqplots II - Checking accurateness of SDM 
====

###1. [Uchida et al](http://pcp.oxfordjournals.org/content/52/4/716.long)
#####Chromosome 1 

- Threshold = 0 

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_1/Rplot.SDM0.png?raw=true)

- Threshold = 2

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_1/Rplot.SDM2.png?raw=true)

- Threshold = 2 and comparison to short expected distribution

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_1/Rplot.SDM2_vs_short.png?raw=true)

#####Chromosome 4

- Threshold = 0 

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/Aw_sup1-2/Variant_calling/sup1_2_4/Rplot.SDM.png?raw=true)

###2. [Galvão et al](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)

- Threshold = 0 
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/OCF2/Rplot.SDM.png?raw=true)

- Threshold = 1
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/OCF2/Rplot.SDM_thres1.png?raw=true)



###3. m mutants 

- Threshold = 0
#####F - Chromosome 1
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/m_mutants/F_chromosome1/Rplot.SDM.png?raw=true)






