
for i in {1..5} 
do
	ruby SNP_distribution_method.rb 1-15Mb/1Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/3Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/3Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/5Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/5Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/7Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/7Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/9Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/9Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/11Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/11Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/13Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/13Mb_A_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/15Mb_$i SDM 0 1
	ruby SNP_distribution_method.rb 1-15Mb/15Mb_A_$i SDM 0 1
done 

for i in {1..5} 
do
	ruby SNP_distribution_method.rb 30Mb/chr1_$i SDM_0 0 1
	ruby SNP_distribution_method.rb 30Mb/chr1_A_$i SDM_0 0 1
	ruby SNP_distribution_method.rb 30Mb/chr1_B_$i SDM_0 0 1

done 

