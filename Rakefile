# desc "make directory"
# task :make do 
# 	sh 'mkdir Sch'
# end 

# desc "Align using bwa"
# task :bwa do
# 	sh 'bwa mem TAIR10_chr_all.fas noadpt6_sch.fq  > Sch/sch.sam'
# end

# desc "Convert sam to bam file"
# task :bam => ["bwa"] do
# 	sh 'samtools view -bS Sch/sch.sam | samtools sort -m 30000000000 - Sch/sch_sorted'
# end

desc "Write pileup file"
task :pileup do
	sh 'samtools mpileup -B -f TAIR10_75/TAIR10_chr_all.fas Parent_BCF2/BCF2_pairs.bam > Parent_BCF2/BCF2_parent.pileup'
end

desc "run VarScan"
task :varscan => ["pileup"] do 
	sh 'java -jar VarScan.v2.3.7.jar mpileup2snp Parent_BCF2/BCF2_parent.pileup --output-vcf 1 > Parent_BCF2/BCF2_parent.vcf'
end
