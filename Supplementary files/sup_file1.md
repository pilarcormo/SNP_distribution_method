Supplementary file 1. Method workflow 
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

1. Index reference sequence 

```
bwa index TAIR10.fa
```

2. Map the reads to reference genome with BWA

```
desc "Align using bwa"
task :bwa  do
      sh 'bwa mem TAIR10.fa paired_R1.fq paired_R2.fq > alignment.sam'
end
```

3. Convert the respective SAM file to BAM file and sort the BAM file using samtools

```
desc "Convert sam to bam file"
task :bam => ["bwa"] do
    sh 'samtools view -bS alignment.sam | samtools sort -m 30000000000 - alignment'
end
```

4. Generate pileup from BAM file

```
desc "Write pileup file"
task :pileup => ["bam"] do
        sh 'samtools mpileup -B -f TAIR10.fa alignment.bam > SNPs.pileup'
end
```

5. Call SNPs using VarScan and record them in a VCF4.1
file 

```
desc "run VarScan"
task :varscan  => ["pileup"] do 
        sh 'java -jar VarScan.v2.3.7.jar mpileup2snp SNPs.pileup --output-vcf 1 
        > SNPs.vcf'
end
```


