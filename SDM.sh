
for i in {1..5} 
do
	echo "$i"
	ruby manage_vcf.rb Galvao $i output_galvao
	ruby manage_vcf.rb Sch $i sch
done 
