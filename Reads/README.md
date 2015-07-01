#SNP calling using NGS reads

We analysed the procedure in existing data sets by examining SNP distribution in recent out-cross (Galvão et al. 2012; Uchida et al. 2014) and back-cross experiments (Allen et al. 2013; Monaghan et al. 2014) in *Arabidopsis thaliana* backgrounds. In all the examples we analysed, homozygous SNPs were normally distributed around the causal mutation. We used the real SNP densities obtained from these experiments to prove the efficiency and accuracy of SDM. The algorithm succeed in the identification of the genomic regions of small size (10-100 kb) containing the causative mutations.


####References

Galvão, Vinicius C; Karl J V Nordström, Christa Lanz, Patric Sulz, Johannes Mathieu, David Posé, Markus Schmid, Detlef Weigel, and Korbinian Schneeberger. 2012. “Synteny-Based Mapping-by-Sequencing Enabled by Targeted Enrichment.” Plant J 71 (3): 517–26. doi:[10.1111/j.1365-313X.2012.04993.x.](http://dx.doi.org/10.1111/j.1365-313X.2012.04993.x)

Uchida, Naoyuki; Tomoaki Sakamoto, Masao Tasaka, and Tetsuya Kurata. 2014. “Identification of EMS-Induced Causal Mutations in Arabidopsis Thaliana by Next-Generation Sequencing.” Methods Mol Biol 1062: 259–70. doi:[10.1007/978-1-62703-580-4_14.](http://dx.doi.org/10.1007/978-1-62703-580-4_14)

Allen, Robert S; Kenlee Nakasugi, Rachel L Doran, Anthony A Millar, and Peter M Waterhouse. 2013. “Facile Mutant Identification via a Single Parental Backcross Method and Application of Whole Genome Sequencing Based Mapping Pipelines.” Front Plant Sci 4: 362. doi:[10.3389/fpls.2013.00362.](http://dx.doi.org/10.3389/fpls.2013.00362)

Monaghan, Jacqueline; Susanne Matschi, Oluwaseyi Shorinola, Hanna Rovenich, Alexandra Matei, Cécile Segonzac, Frederikke Gro Malinovsky, et al. 2014. “The Calcium-Dependent Protein Kinase CPK28 Buffers Plant Immunity and Regulates BIK1 Turnover.” Cell Host Microbe 16 (5): 605–15. doi:[10.1016/j.chom.2014.10.007.](http://dx.doi.org/10.1016/j.chom.2014.10.007)

Method workflow 
===

####Quality filtering for paired-end reads

```
FastQC/fastqc reads_R1.fq
FastQC/fastqc reads_R2.fq
```

```
java -jar Trimmomatic-0.33/trimmomatic-0.33.jar PE reads_R1.fq reads_R2.fq paired_R1.fq unpaired_R1.fq paired_R2.fq unpaired_R2.fq  SLIDINGWINDOW:4:20 MINLEN:70"
```

####Command line parameters used for paired-end read mapping and SNP calling

The reference genome was indexed before running the alingment and SNP calling pipeline. 

```
bwa index TAIR10.fa
```

```
desc "Align using bwa"
task :bwa  do
      sh 'bwa mem TAIR10.fa paired_R1.fq paired_R2.fq > alignment.sam'
end
```
```
desc "Convert sam to bam file"
task :bam => ["bwa"] do
    sh 'samtools view -bS alignment.sam | samtools sort -m 30000000000 - alignment'
end
```
```
desc "Write pileup file"
task :pileup => ["bam"] do
        sh 'samtools mpileup -B -f TAIR10.fa alignment.bam > SNPs.pileup'
end
```
```
desc "run VarScan"
task :varscan  => ["pileup"] do 
        sh 'java -jar VarScan.v2.3.7.jar mpileup2snp SNPs.pileup --output-vcf 1 > SNPs.vcf'
end
```

VCF files obtained are in [https://github.com/pilarcormo/SNP_distribution_method/tree/master/Reads](https://github.com/pilarcormo/SNP_distribution_method/tree/master/Reads)

###Filtering background SNPs 


Run [manage_vcf.rb](https://github.com/pilarcormo/SNP_distribution_method/blob/master/manage_vcf.rb)

 ``ruby manage_vcf.rb (1) (2) (3) (4) -(5)-``

1. location folder
2. Variants in VCF 4.1 file
3. the chromosome we want to analyse (1, 2, 3, 4, 5) 
4.  a method for managing the vcf file: use either 'cutting_vcf' or 'filter_vcf'. 
If the second option is used we also need 
5. a parental (or long distant ecotype) variants in VCF 4.1 file. 

The **cutting_vcf** option will create in individual VCF file for each chromosome. The **filter_vcf** option removes the SNPs from the mutant VCF file if they are also present in the parental VCF provided.

Also, it simplifies the output VDF. In the INFO field each heterozygous SNP will have AF = 0.5, while homozygous SNPs will have AF = 1.0 to facilitate the modelling of SDM. It will also create text files for homozygouys and heterozygous SNP positions.


###Removing the centromeres 
Run [remove_cent.rb](https://github.com/pilarcormo/SNP_distribution_method/blob/master/remove_cent.rb)

 ```ruby remove_cent.rb chromosome folder_containing_SNPs_files``` 
 
It takes the SNP lists obtained from the VCF file and remove the SNP positions that are caused by the high variability in the centromeric regions. For now, it is defined for *Arabidopsis thaliana* only.

```
centromere = {"chr1" => [15086545-3950000, 15086545+3950000],"chr2" => [3608429- 1500000, 3608429+ 1500000], "chr3" => [14209452- 1500000, 14209452+1500000], "chr4" => [3956521- 1400000, 3956521+1400000], "chr5" => [11725524-500000, 11725524+500000]}

```

###Creating model genomes based on real SNP densities from reverse genetic screens 

To create model genomes with real SNP densities obtained from variant calling with NGS reads, we can run [model_genome_real_hpc.rb](https://github.com/pilarcormo/SNP_distribution_method/model_genome_real_hpc.rb)

```ruby model_genome_real_hpc.rb dataset_name contig_size folder_containing_SNPs_files chromosome```


Example:
```
ruby model_genome_real_hpc.rb B_nocen_chr5_10kb 10000  /Users/morenop/SNP_distribution_method/Reads/m_mutants/B_chromosome5/interesting_5/ 5
```

Instead of using an ideal normal distribution for the homozygous SNPs and a uniform distribution for the heterozygous SNPs, the SNPs present in the VCF files after the SNP calling are used as input. In each case, the chromosome where the mutation is described is divided in contigs of size expecified.

The model genomes created before removing the centromere are available at 
[https://github.com/pilarcormo/SNP_distribution_method/tree/master/arabidopsis_datasets/No_centromere](https://github.com/pilarcormo/SNP_distribution_method/tree/master/arabidopsis_datasets/No_centromere).  The "_1kb", "_10kb", etc in the name indicates the minimum contig size used (contig size will oscilate between this value and its double).

The model genomes created after removing the centromere are available at [https://github.com/pilarcormo/SNP_distribution_method/tree/master/arabidopsis_datasets/Centromere](https://github.com/pilarcormo/SNP_distribution_method/tree/master/arabidopsis_datasets/Centromere).

##Probability plots 
To check if the homozygous SNP distributions obtained from SNP calling in back-cross and out-cross experiments correlate to a normal distribution, we used QQ-plots. Results and R code available at [https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/qqplot.md](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/qqplot.md)

Then, the correlation (r2), standard deviation, kurtosis and skewness of the homozygous SNP density was measured. Results and R code available at [https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/qqplot.md](https://github.com/pilarcormo/SNP_distribution_method/blob/master/Reads/kurtosis.md)


###Project dependencies

1. Ruby >= 2

2. Ruby gems:

	- bio >= 1.4.3.0001
	- bio-samtools >= 2.2.0
	- rinruby >= 2.0.3

3. R >= 3.1.1


