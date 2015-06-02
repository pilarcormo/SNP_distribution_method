SNP Distribution Method
======

Whole genome sequencing using Next-Generation Sequencing (NGS) technologies offers a unique opportunity to study genetic variations. However, mapping the mutations responsible for phenotypes is generally a tedious and time-consuming process. In the last few years, researchers have developed user-friendly tools to identify mutations, yet they are not applicable to organisms with non-sequenced genomes. 
We aim to develop a new tool that can locate causative SNP mutations in backcrossing experiments without a reference genome. Forward genetic screens offer the opportunity to experimentally induce a phenotype of interest in a target organism. 
For our purpose, a mutant individual is crossed to a non-mutant line followed by one generation of backcrossing of the heterozygous F1 samples. The offspring of this second cross gives rise to a recombinant population which will segregate for the mutant phenotype. Due to linkage to the phenotype altering SNP, the remaining linked homozygous SNPs will be distributed around it generating a high homozygous SNP density in this non-recombinant genomic region. Consequently, the homozygous/heterozygous SNP ratio will be higher in the area where the SNP of interest is located.
To emulate real data, I created a model genome with the expected SNP density. Then, I split it into fragments that will imitate the contigs assembled from NGS reads. My tool will need to rearrange these contigs based on their SNP density in order to find the causative mutation. The method I have developed, called SNP Distibution Method (SDM), first groups the contigs by their normalised SNP density and then it arranges them following the ideal distribution. As it was expected, we had to eliminate the contigs with a low homozygous SNP density from the analysis since they distorted the accuracy of the approach. This method was able to effectively identify the artificially created causative mutation when using the chromosome 1 from A. thaliana as model genome. 

######How does it work?

If we have have a fasta file with the unordered contigs and a VCF file with the number of hm and ht SNPs per contig and we know the SNPs are clustered together around the causative mutation, then we can order the contigs in the first file by taking the 2 lowest values of SNP frecuency and putting them at both ends of the distribution. By repeting this step, we will potentially get to the highest value = causative mutation = peak of the distribution. As contigs do not have a constant length, for the method to be effective, instead of using the absolute number of SNPs, I use a relative number by dividing the SNP density per contig by the number of nucleotides (length).

###Creating a model genome

Running ```ruby create_model_genome.rb dataset_name genome_length contig_size``` will generate a new model genome in SNP_distribution_method/arabidopsis_datasets/dataset_name. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0. The variables hm_r and ht_r contain the R code needed to create the model homozygous and heterozygous SNP distributions respectively. The variable contig_size provides the minimum size for contigs, where the maximum size is double this value, and each contig's size is randomly chosen within this range. Information obtained from [https://github.com/edwardchalstrey1/fragmented_genome_with_snps](https://github.com/edwardchalstrey1/fragmented_genome_with_snps)

To create model genomes with real SNP densities obtaind from variant calling with NGS reads, we can run ```ruby model_genome_real_hpc.rb dataset_name contig_size folder_containing_SNPs_files chromosome```


###Filtering background SNPs 

Run ``ruby manage_vcf.rb (1) (2) (3) (4) -(5)-``

1. location folder
2. Variants in VCF 4.1 file
3. the chromosome we want to analyse (1, 2, 3, 4, 5) 
4.  a method for managing the vcf file: use either 'cutting_vcf' or 'filter_vcf'. 
If the second option is used we also need 
5. a parental (or long distant ecotype) variants in VCF 4.1 file. 

The **cutting_vcf** option will create in individual VCF file for each chromosome. The **filter_vcf** option removes the SNPs from the output VCF file if they are also present in the parental VCF provided.

The output VCF files will be more simple.  In the INFO field each heterozygous SNP will have AF = 0.5, and homozygous will have AF = 1.0
to facilitate the modelling of SDM. It will also create text files for homozygouys and heterozygous SNP positions.


###Removing the centromeres 
Run ```ruby remove_cent.rb chromosome folder_containing_SNPs_files```. It takes the SNP lists and remove the SNP positions that are caused by the high variability in the centromeric regions. For now, it is defined for *Arabidopsis thaliana* only.

###Runing SDM

Run ```ruby SNP_distribution_method_variation.rb (1) (2) (3) (4) (5)```

1. dataset_name
2. name for the output folder
3.  Integer 1/0 = filtering step on/off. If it's on, contigs which a ratio below the threshold will be discarded.
4.  Float (1, 0.1, 0.01...). Factor to calculate the ratio.
5.   back/out (kind of cross)

The dataset folder is obtained by running the model_genome explained above. It will generate a FASTA file with the correctly ordered fragments, another FASTA file the shuffled fragments, text files with the list of homozygous and heterozygus SNPs and a VCF file with the SNPs. 

The output after running SDM will be a new FASTA file with the suggested order of contigs, and inside the output folder we will have: 

- Text files for homozygous and heterozygous SNPs after sorting step (perm_hm and perm_ht) and for the hypothetical ratio (hyp_ratio). 
- Text files for homozygous and heterozygous SNPs after pre-filtering step (hm_snps_short, ht_snps_short) in the correctly ordered genome and  for the ratios in the correctly ordered fragments (ratio).
- A plot for the hypothetical SNP densities and ratio
- A QQ-plot comparing the correlation of the obtained SNP density and the expected normal distribution
- A plot comparing the real ratio distribution and the ratio distribution after running SDM 
- A text file (mutation) with the following information:
	 - Length of the group of contigs that form the peak of the distribution
	 - Contigs and positions where the mutation is likely to be found 	 

###R scripts 

All the density plots, QQplots and histograms are described in the R scripts at [https://github.com/pilarcormo/SNP_distribution_method/R_scripts](https://github.com/pilarcormo/SNP_distribution_method/R_scripts)

###Project dependencies

1. Ruby >= 2.0.0

2. Ruby gems:

	- bio >= 1.4.3.0001
	- bio-samtools >= 2.2.0
	- rinruby >= 2.0.3

3. R >= 3.1.1


