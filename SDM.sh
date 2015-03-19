
for i in {0..5} 
do
	echo "$i"
	ruby SNP_distribution_method_variation.rb BCF2_4_v2 Perm_1903_ratio_$i $i
done 
