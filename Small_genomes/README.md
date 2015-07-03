
#Modelling SDM using small genomes


###Creating a model genome

For our purpose, a mutant individual should be crossed to a non-mutant parenta line (either the same ecotype -backcross- or a distant ecotype -outcross-). The offspring of this cross would give rise to a recombinant population which will segregate for the mutant phenotype. Due to linkage to the phenotype altering SNP, the remaining linked homozygous SNPs will be distributed around it generating a high homozygous SNP density in this non-recombinant genomic region. Consequently, the homozygous/heterozygous SNP ratio will be higher in the area where the SNP of interest is located.

To identify the causal mutation in the model genomes, we used idealised SNP distributions. We observed that Homozygous SNPs frequency around the causative mutation followed a normal distribution, with the causal mutation  in the middle of the distribution.  Heterozygous SNPs followed a uniform distribution, with no specific contribution to the general SNP density.

Running ```ruby create_model_genome.rb dataset_name genome_length contig_size``` will generate a new model genome in SNP_distribution_method/arabidopsis_datasets/{dataset_name}. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0. The variables hm_r and ht_r contain the R code needed to create the model homozygous and heterozygous SNP distributions respectively. The variable contig_size provides the minimum size for contigs, where the maximum size is double this value, and each contig's size is randomly chosen within this range. Information obtained from [https://github.com/edwardchalstrey1/fragmented_genome_with_snps](https://github.com/edwardchalstrey1/fragmented_genome_with_snps). 

To create these model genomes (1-15Mb), a SNP density of 1 SNP/500 bp was used using the following R code:

```
# Create the lists of homozygous and heterozygous SNPs
hm_r = "hm <- rnorm(#{snp}, #{size/2}, #{snp*2})" 
ht_r = "ht <- runif(#{snp}, 1, #{size})"   
```
where size is the genome size in bp and snp is ```snp = (genome_size/1000)*2```

The folders without letter contain genomes divided in 700 contigs and the folders with an "A" contain genomes with 1300 contigs. 5 replicates were created for each genome size and contig size (1-5). They can be found at [https://github.com/pilarcormo?Small_genomes/arabidopsis_datasets/1-15Mb](https://github.com/pilarcormo/Small_genomes/arabidopsis_datasets/1-15Mb). The info.txt file contains the homozygous and heterozygous SNPs densities used to create the genome and the contig size used in each case. 

Then, I also created [https://github.com/pilarcormo?Small_genomes/arabidopsis_datasets/30Mb](https://github.com/pilarcormo/Small_genomes/arabidopsis_datasets/30Mb)


###Causal mutation not located in the mean
The 30 Mb genomes under the name _C and _E were created with the homozygous SNP mutation moved to the right tail and the left tail (respectively). Instead of defining the mean of the rnorm distribution in the middle of the distribution, this value was shifted around a 20% to the right and left.

```
# Create the lists of homozygous and heterozygous SNPs
hm_r = "hm <- rnorm(#{snp}, mean, #{snp*2})" # Causative SNP at/near 10000
ht_r = "ht <- runif(#{snp}, 1, #{size})"   # Genome length of 10000
```
i = 1..5
<table>
 <tr><th>Name <th>Contig length</th> <th>Mean of rnorm distribution</th>
 <tr><th>arabidopsis_datastets/30Mb/chr1_C_i <th>10000 </th> <th>20000000</th>
 <tr><th>arabidopsis_datastets/30Mb/chr1_E_i <th>10000 </th> <th>10000000</th>
</table>

A [pre-filtering step](https://github.com/pilarcormo/Small_genomes/arabidopsis_datasets/Analyse_effect_ratio) based on the homozygous to heterozygous SNPs ratio was incorporated to SDM to remove the contigs that are surrounding the causal mutation (those located in the tails of the normal distribution of homozogyoys SNPs). This filtering step removes the contigs with a hom/het ratio below a given percentage of the maximum ratio in the assembly.


###Runing SDM


Look at SDM.sh for the shell script used to run SDM on the model genomes with different sizes and contig lenghts. 

Run ```ruby SNP_distribution_method_variation.rb (1) (2) (3)```

1. **Input dataset folder** containing the input files -Contigs in FASTA file and SNPs in a VCF file-
2. name for the **output folder**
3. **Threshold = provide one of the following 0, 1, 5, 10, 20.**
	- 0 -> filtering step off. 
	- larger than 0 -> filtering step on.  If the filtering step is required, the threshold astringency is provided as an integer (1, 5, 10, 20). Each integer represents the percentage of the maximum ratio below which a contig will be discarded. In example, if 1 is specified, SDM will discard those contigs with a ratio falling below 1% of the maximum ratio while a value of 20 is more astringent  will discard those contigs with a ratio falling below 20% of the maximum ratio. 

To test SDM, the input dataset folder (1) is obtained by running the model_genome explained above. It will generate a FASTA file with the correctly ordered fragments, another FASTA file the shuffled fragments, text files with the list of homozygous and heterozygus SNPs and a VCF file with the SNPs. 

The output after running SDM will be a new FASTA file with the suggested order of contigs, and inside the output folder we will have. The model genomes generated to test SDM are in 

######Deviation from causal mutation

```
library(ggplot2)
library(grid)
library(gridExtra)
```
```
deviations <- read.csv("~/SNP_Distribution_method/Small_genomes/arabidopsis_datasets/1-15Mb.csv")
deviations_30 <- read.csv("~/SNP_Distribution_method/Small_genomes/arabidopsis_datasets/30Mb.csv")
deviations$Genome_size <- factor(deviations$Genome_size, labels = c("1", "3", "5", "7", "9", "11", "13", "15"))
deviations$Contigs <- factor(deviations$Contigs, labels = c("1300", "700"))
deviations_30$Genome_size <- factor(deviations_30$Genome_size, labels = c("30"))
deviations_30$Contigs <- factor(deviations_30$Contigs, labels = c("4000", "2000", "1000"))
```
```
Palette <- c('brown3',"royalblue2")
Palette2 <- c('green4',"green3", "greenyellow")
g <- ggplot(deviations, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) + geom_point(width = .2) + ylim(0, 3.2) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette) + theme_bw()
h <- ggplot(deviations_30, aes(x = Genome_size, y = Deviation, shape = Contigs, colour = Contigs)) +  geom_point(width = .2) + ylim(0, 3.2) + labs(x = "Genome size (Mb)",  y = "% of deviation") + theme_bw() + scale_colour_manual(values=Palette2) + theme_bw()
gh <- grid.arrange(g , h, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)
```

![Image](Rplot.jitter_order.png)

###Project dependencies

1. Ruby >= 2.0.0

2. Ruby gems:

	- bio >= 1.4.3.0001
	- bio-samtools >= 2.2.0
	- rinruby >= 2.0.3

3. R >= 3.1.1


