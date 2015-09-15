source ruby-2.0.0
source R-3.1.0


# bsub -q TSL-Test128 -o BCF2_real.txt "ruby model_genome_new.rb BCF2_1 10000 1"
# bsub -q TSL-Test128 -o BCF2_real.txt "ruby model_genome_new.rb BCF2_2_v2 10000 2"
# bsub -q TSL-Test128 -o BCF2_real.txt "ruby model_genome_new.rb BCF2_3 10000 3"
# bsub -q TSL-Test128 -o int2oc_real.txt "ruby model_genome_new.rb OC_chr2_1kb 1000 2 OF/Interesting_2"
# bsub -q TSL-Test128 -o int2oc_real.txt "ruby model_genome_new.rb OC_chr2_10kb 10000 2 OF/Interesting_2"
# bsub -q TSL-Test128 -o int2oc_real.txt "ruby model_genome_new.rb OC_chr2_100kb 100000 2 OF/Interesting_2"

# bsub -q TSL-Test128 -o mob.txt "ruby model_genome_new.rb mob1_chr5_1kb 1000 5 mob62_interesting_5"
# bsub -q TSL-Test128 -o mob.txt "ruby model_genome_new.rb mob1_chr5_10kb 10000 5 mob62_interesting_5"
# bsub -q TSL-Test128 -o mob.txt "ruby model_genome_new.rb mob1_chr5_100kb 100000 5 mob62_interesting_5"

# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb bcf2_nocen_chr3_100kb 100000  new_snps/bcf2_3"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb bcf2_last_chr3_100kb 100000 new_snps/bcf2_3 3"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb ocf2_lat_chr2_100kb 100000 new_snps/ocf2 2"



# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model2_sup1_nocen_ch2_10kb 2000 new_snps/filter2_chromosome4 4"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model5_sup1_nocen_ch5_10kb 5000 new_snps/filter2_chromosome4 4"

# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model2_bcf2_nocen_ch2_10kb 2000 new_snps/bcf2_3 3"
# # bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model5_bcf2_nocen_ch5_10kb 5000 new_snps/bcf2_3 3"

# # bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model2_ocf2_nocen_ch2_10kb 2000 new_snps/ocf2 2"
# # bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model5_ocf2_nocen_ch5_10kb 5000 new_snps/ocf2 2"

# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model2_C_nocen_chr5_2kb 2000 new_snps/C_5 5"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model5_C_nocen_chr5_5kb 5000 new_snps/C_5 5"

bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model2_B_nocen_chr5_2kb 2000 new_snps/B_5 5"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb model5_B_nocen_chr5_5kb 5000 new_snps/B_5 5"




# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb C_nocen_chr5_10kb 10000 new_snps/C_5 5"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb ocf2_nocen_chr2_10kb 10000 new_snps/ocf2 2"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb bcf2_nocen_chr3_10kb 10000 new_snps/bcf2_3 3"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb sup1_nocen_chr4_10kb 10000 new_snps/sup1_4 4"



# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb BC_nocen_chr5_100kb 100000 new_snps/BC 5"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb BC_nocen_chr5_100kb 100000 new_snps/BC 5"



# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb E_nocen_whole_100kb 100000 new_snps/C_wholegenome"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb F_nocen_whole_100kb 100000 new_snps/F_wholegenome"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb ocf2_nocen_whole_100kb 100000 new_snps/ocf2_wholegenome"
# bsub -q TSL-Test128 -o nocen.txt "ruby model_genome_new.rb sup1_nocen_whole_100kb 100000 new_snps/sup1_wholegenome"



# bsub -q TSL-Test128 -o int2oc_real.txt "ruby model_genome_new.rb OC_chr2 10000 2 OF/Interesting_2"
# bsub -q TSL-Test128 -o BCF2_real.txt "ruby model_genome_new.rb BCF2_5 10000 5"

# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_3Mb Perm2"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_5Mb Perm"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_7Mb Perm"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_9Mb Perm"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_11Mb Perm"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_13Mb Perm"
# bsub -q TSL-Test128 -o G.txt "ruby find_causal_mutation.rb G_15Mb Perm"


# bsub -q TSL-Prod128 -o chr1.txt "ruby SNP_distribution_method.rb Chromosome1 Perm"

# bsub -q TSL-Prod128 -o chr1.txt "ruby find_causal_mutation.rb Chromosome1 Perm"


# for i in {1..5} 
# do
# 	bsub -q TSL-Test128 -o 1ch_1902.txt "ruby SNP_distribution_method.rb chr1_E_$i Perm_1902"
# 	# bsub -q TSL-Test128 -o 1ch_1902.txt "ruby SNP_distribution_method.rb chr1_C_$i Perm_1902"
# 	# bsub -q TSL-Test128 -o 1ch_2601.txt "ruby SNP_distribution_method.rb chr1_C_$i Perm"
# 	# bsub -q TSL-Prod128 -o 1ch_2601.txt "ruby SNP_distribution_method.rb chr1_D_$i Perm"
# 	# bsub -q TSL-Prod128 -o 1ch_2601.txt "ruby SNP_distribution_method.rb chr1_E_$i Perm"
# done 



# for i in {1..5} 
# do
# # 	bsub -q TSL-Prod128 -o 1ch.txt "ruby SNP_distribution_method.rb chr1_$i Perm"
# # 	bsub -q TSL-Test128 -o 1ch.txt "ruby SNP_distribution_method.rb chr1_A_$i Perm"
# 	bsub -q TSL-Prod128 -o 1ch.txt "ruby SNP_distribution_method.rb chr1_B_$i Perm"
# done 



# for i in {1..5} 
# do
# 	# bsub -q TSL-Prod128 -o 1ch.txt "ruby find_causal_mutation.rb chr1_$i Perm"
# # 	bsub -q TSL-Test128 -o 1ch.txt "ruby find_causal_mutation.rb chr1_A_$i Perm"
# 	bsub -q TSL-Test128 -o 1ch.txt "ruby find_causal_mutation.rb chr1_B_$i Perm"
# done 