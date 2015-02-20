SNP Distribution Method
======

Whole genome sequencing using Next-Generation Sequencing (NGS) technologies offers a unique opportunity to study genetic variations. However, mapping the mutations responsible for phenotypes is generally a tedious and time-consuming process. In the last few years, researchers have developed user-friendly tools to identify mutations, yet they are not applicable to organisms with non-sequenced genomes. 
We aim to develop a new tool that can locate causative SNP mutations in backcrossing experiments without a reference genome. Forward genetic screens offer the opportunity to experimentally induce a phenotype of interest in a target organism. 
For our purpose, a mutant individual is crossed to a non-mutant line followed by one generation of backcrossing of the heterozygous F1 samples. The offspring of this second cross gives rise to a recombinant population which will segregate for the mutant phenotype. Due to linkage to the phenotype altering SNP, the remaining linked homozygous SNPs will be distributed around it generating a high homozygous SNP density in this non-recombinant genomic region. Consequently, the homozygous/heterozygous SNP ratio will be higher in the area where the SNP of interest is located.
To emulate real data, I created a model genome with the expected SNP density. Then, I split it into fragments that will imitate the contigs assembled from NGS reads. My tool will need to rearrange these contigs based on their SNP density in order to find the causative mutation. The method I have developed, called SNP Distibution Method (SDM), first groups the contigs by their normalised SNP density and then it arranges them following the ideal distribution. As it was expected, we had to eliminate the contigs with a low homozygous SNP density from the analysis since they distorted the accuracy of the approach. This method was able to effectively identify the artificially created causative mutation when using the chromosome 1 from A. thaliana as model genome. 


####How does it work?

If we have have a fasta file with the unordered contigs and a VCF file with the number of hm and ht SNPs per contig and we know the SNPs are clustered together around the causative mutation, then we can order the contigs in the first file by taking the 2 lowest values of SNP frecuency and putting them at both ends of the distribution. By repeting this step, we will potentially get to the highest value = causative mutation = peak of the distribution. As contigs do not have a constant length, for the method to be effective, instead of using the absolute number of SNPs, I use a relative number by dividing the SNP density per contig by the number of nucleotides (length).


####What's next?
We now aim to use real NGS data to check the usefulness of SDM in a real-life experiment. Then, we will need to develop a tool that is fast and accessible in non-bionformatics environments. 