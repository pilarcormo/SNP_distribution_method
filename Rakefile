desc "make directory"
task :make do 
	sh 'mkdir Sch'
end 

desc "Align using bwa"
task :bwa do
	sh 'bwa mem TAIR10_chr_all.fas noadpt6_sch.fq  > Sch/sch.sam'
end

desc "Convert sam to bam file"
task :bam => ["bwa"] do
	sh 'samtools view -bS Sch/sch.sam | samtools sort -m 30000000000 - Sch/sch_sorted'
end

desc "Write pileup file"
task :pileup => ["bam"] do
	sh 'samtools mpileup -B -f TAIR10_chr_all.fas Sch/sch_sorted.bam > Sch/sch.pileup'
end

desc "run VarScan"
task :varscan => ["pileup"] do 
	sh 'java -jar VarScan.v2.3.7.jar mpileup2snp Sch/sch.pileup --output-vcf 1 > Sch/sch.vcf'
end
