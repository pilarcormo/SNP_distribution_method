
1. Building individual VCF files for each chromosome

ruby manage_vcf.rb Reads/Aw_sup1-2/Variant_calling sup1 sup1_2_$i $i cutting_vcf
ruby manage_vcf.rb Reads/Aw_sup1-2/Variant_calling sup1 sup1_2_$i $i cutting_vcf

for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2 BCF2 BCF2_chromosome$i $i cutting_vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2/BCF2_chromosome$i $i interesting_$i $i filter_vcf reads/BCF2/BCF2_parent/BCF2_parent.vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/m_mutants B B_chromosome$i $i cutting_vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/m_mutants/B_chromosome$i $i interesting_$i $i filter_vcf reads/m_mutants/parent/parent.vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/m_mutants C C_chromosome$i $i cutting_vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/m_mutants/C_chromosome$i $i interesting_$i $i filter_vcf reads/m_mutants/parent/parent.vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2 BCF2 BCF2_chromosome$i $i cutting_vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2/BCF2_chromosome$i $i interesting_$i $i filter_vcf reads/BCF2/BCF2_parent/BCF2_parent.vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2 BCF2 BCF2_chromosome$i $i cutting_vcf
done


for i in {1..5} 
do
ruby manage_vcf.rb Reads/BCF2_2/BCF2_chromosome$i $i interesting_$i $i filter_vcf reads/BCF2/BCF2_parent/BCF2_parent.vcf
done


Reads/BCF2 BCF2 BCF2_chromosome$i $i cutting_vcf
filter_vicf

ruby manage_vcf.rb reads/Alpina arabis_bc1f2 3 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc1f2 4 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc1f2 5 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc1f2 6 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc1f2 7 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc1f2 8 cutting_vcf

ruby manage_vcf.rb reads/Alpina arabis_bc2f2 1 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 2 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 3 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 4 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 5 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 6 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 7 cutting_vcf
ruby manage_vcf.rb reads/Alpina arabis_bc2f2 8 cutting_vcf