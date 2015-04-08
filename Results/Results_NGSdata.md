
Results
====

###[BCF2 reads](http://bioinfo.mpipz.mpg.de/shoremap/SHOREmap_v3.0.html)

In here I'm showing the overlapping plots with both densities. In red, the hm SNP density and in blue, the ht SNP density. By doing this, I can start guessing how the hm/ht ratio distribution will look like. 

#####Chromosome 1 
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome1/Rplot.hmhtdensities_magblue.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome1/Rplot.hist%26den.png?raw=true)


#####Chromosome 2 
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome2/Rplot.hmhtdensities_magblu.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome2/Rplot.hist%26den.png?raw=true)

#####Chromosome 3
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome3/Rplot.hmht.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome3/Rplot.hist%26den.png?raw=true)


#####Chromosome 4
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome4/Rplot.hmhtdensities_magblu.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome4/Rplot.hist%26den.png?raw=true)
#####Chromosome 5 
![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome5/Rplot.hmhtdensities_magblu.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome5/Rplot.hist%26den.png?raw=true)

###[Galvao/OCF2 reads](http://bioinfo.mpipz.mpg.de/shoremap/SHOREmap_v3.0.html)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome1/Rplot.hmht.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome2/Rplot.hmht.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome3/Rplot.hmht.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome4/Rplot.hmht.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome5/Rplot.hmht.png?raw=true)



###Causative SNP

The results for BCF2 are explained in [this paper by Allen et al](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772335/#SM3). They found the causative SNP at the beginning of chromosome 3  (microRNA pathway mutant). However, with my SNP data I don't see any high SNP density on this area of the genome. I plot the normalised ratio against the genomic position:

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/BCF2/BCF2_chromosome3/Rplot.ratio.png?raw=true)

For [OCF2](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9), they found the mutation responsible for the late flowering in chromosome 2. However, with my vcf file I cannnot see any interesting SNP distribution around this area. In the paper they also say that **"along chromosome 2, the increases or decreases in AFEs around the causal mutation were not monotonic, as expected in mapping populations for a single Mendelian trait. This suggests that there is a considerable difference between real allele frequencies and AFEs."** They calculated the Allele Frequency estimators (AFEs) that are the percentage of reads supporting one of the parental alleles. I plot the ratios observed again:


![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome2/Rplot.ratio.png?raw=true)

![Image](https://github.com/pilarcormo/SNP_distribution_method/blob/master/OF/chromosome2/Rplot.10-20Mb.png?raw=true)


