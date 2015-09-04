####mob mutants
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
####OCF2
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

####BCF2

for i in {1..5} 
do
	ruby manage_vcf.rb Reads/BCF2 BCF2.vcf BCF2_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 BCF2_parent.vcf chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/BCF2 BCF2_chromosome$i/chromosome$i.vcf interesting_$i $i filter_vcf BCF2_parent/chromosome$i/chromosome$i.vcf
done 

####sup1

for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 sup1.vcf sup1_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 colT.vcf Parental/colT_chromosome$i $i cutting_vcf
	ruby manage_vcf.rb Reads/Aw_sup1-2 WsT.vcf Parental/Ws_chromosome$i $i cutting_vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 sup1_chromosome$i/chromosome$i.vcf filter1_chromosome$i $i filter_vcf Parental/Ws_chromosome$i/chromosome$i.vcf
done 
for i in {1..5} 
do
	ruby manage_vcf.rb Reads/Aw_sup1-2 filter1_chromosome$i/chromosome$i.vcf filter2_chromosome$i $i filter_vcf Parental/ColT_chromosome$i/chromosome$i.vcf
done 