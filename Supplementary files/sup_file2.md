
Supplementary file 2. Filtering workflow 
===
###Parental filtering

####mob mutants
```
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/m_mutants B.vcf B_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/m_mutants C.vcf C_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/m_mutants/parent parent.vcf chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/m_mutants B_chromosome$i/chromosome$i.vcf interesting_$i $i filter_vcf Reads/m_mutants/parent/chromosome$i/chromosome$i.vcf
	ruby manage_vcf.rb Reads/m_mutants C_chromosome$i/chromosome$i.vcf interesting_$i $i filter_vcf Reads/m_mutants/parent/chromosome$i/chromosome$i.vcf
done 
```
####OCF2

```
unzip Reads/OCF2/OF_output25vcf.zip
unzip Reads/OCF2/Ler/OC_parent.vcf.zip
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/BCF2 OF/OF_output25.vcf OCF2_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/OCF2/Ler OCF2/OCF2_parent/OC_output.vcf chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/OCF2 OCF2_chromosome$i/chromosome$i.vcf interesting_$i $i filter_vcf Ler/chromosome$i/chromosome$i.vcf
done 
```

####BCF2

```
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/BCF2 BCF2.vcf BCF2_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 BCF2_parent.vcf chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/BCF2 BCF2_chromosome$i/chromosome$i.vcf interesting_$i $i filter_vcf BCF2_parent/chromosome$i/chromosome$i.vcf
done 
```

####sup1

```
unzip Reads/Aw_sup1-2/vcfs.zip
unzip Reads/Parental/vcfs.zip
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 vcfs/sup1.vcf sup1_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 Parental/vcfs/colT.vcf Parental/colT_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 Parental/vcfs/WsT.vcf Parental/Ws_chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 sup1_chromosome$i/chromosome$i.vcf 	filter1_chromosome$i $i filter_vcf Parental/Ws_chromosome$i/chromosome$i.vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 filter1_chromosome$i/chromosome$i.vcf filter2_chromosome$i $i filter_vcf Parental/ColT_chromosome$i/chromosome$i.vcf
done 
```

###Centromere filtering
```
for i in {1..5} 
do
	ruby remove_cent.rb $i Aw_sup1-2/filter2_chromosome$i
	ruby remove_cent.rb $i BCF2/BCF2_chromosome$i
	ruby remove_cent.rb $i OCF2/OCF2_chromosome$i
	ruby remove_cent.rb $i B/B_chromosome$i
	ruby remove_cent.rb $i C/C_chromosome$i
done
```
###SNP density analysis
```
for i in {1..5}; do ruby snp_density.rb $i BCF2	BCF2_chromosome$i interesting_$i BCF2; done 
for i in {1..5}; do ruby snp_density.rb $i OCF2	OCF2_chromosome$i Interesting_$i OCF2; done 
for i in {1..5}; do ruby snp_density.rb $i m_mutants B_chromosome$i interesting_$i mob1; done 
for i in {1..5}; do ruby snp_density.rb $i m_mutants C_chromosome$i interesting_$i mob2; done 
for i in {1..5}; do ruby snp_density.rb $i Aw_sup1-2 sup1_chromosome$i ../filter2_chromosome$i sup1; done 
```

###Model genomes with real densities

```
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/sup1_nocen_chr4_2kb	2000 arabidopsis_datasets/SNP_densities/sup1_4 4
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/sup1_nocen_chr4_5kb	5000 arabidopsis_datasets/SNP_densities/sup1_4 4
ruby model_genome_real_hpc.rb No_centromere/10kb_contig/sup1_nocen_chr4_10kb	10000 arabidopsis_datasets/SNP_densities/sup1_4 4
```
```
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/bcf2_nocen_ch2_2kb	2000 arabidopsis_datasets/SNP_densities/bcf2_3 3
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/bcf2_nocen_ch5_5kb	5000 arabidopsis_datasets/SNP_densities/bcf2_3 3
ruby model_genome_real_hpc.rb No_centromere/10kb_contig/bcf2_nocen_ch5_10kb	10000 arabidopsis_datasets/SNP_densities/bcf2_3 3
```
```
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/ocf2_nocen_chr3_2kb	2000 arabidopsis_datasets/SNP_densities/ocf2_2 2
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/ocf2_nocen_chr3_5kb	5000 arabidopsis_datasets/SNP_densities/ocf2_2 2
ruby model_genome_real_hpc.rb No_centromere/10kb_contig/ocf2_nocen_chr3_10kb	10000 arabidopsis_datasets/SNP_densities/ocf2_2 2
```
```
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/C_nocen_chr5_2kb	2000 arabidopsis_datasets/SNP_densities/C_5 5
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/C_nocen_chr5_5kb	5000 arabidopsis_datasets/SNP_densities/C_5 5
ruby model_genome_real_hpc.rb No_centromere/10kb_contig/C_nocen_chr5_10kb	10000 arabidopsis_datasets/SNP_densities/C_5 5
```
```
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/B_nocen_chr5_2kb	2000 arabidopsis_datasets/SNP_densities/B_5 5
ruby model_genome_real_hpc.rb No_centromere/2-5kb_contig/B_nocen_chr5_5kb	5000 arabidopsis_datasets/SNP_densities/B_5 5
ruby model_genome_real_hpc.rb No_centromere/10kb_contig/B_nocen_chr5_10kb	10000 arabidopsis_datasets/SNP_densities/B_5 5
```