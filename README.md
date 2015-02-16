SNP Distribution Method
======



If we have have a file with the unordered contigs and a file with the number of hm and ht SNPs per contig and we know the SNPs are cluster around the causative mutation, we can order the contigs in the first file by taking the 2 lowest values of SNP frecuency and putting them at both ends of the distribution. By repeting this step, we will potentially get to the highest value = causative mutation = peak of the distribution. As contigs do not have a constant length, for the method to be effective, instead of using the absolute number of SNPs, I use a relative number by dividing the SNP density per contig by the number of nucleotides (length).
