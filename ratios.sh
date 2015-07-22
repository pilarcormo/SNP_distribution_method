
#####OCF2
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb SDM_0 0 1 out ###filtering step off
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb SDM_1 1 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb SDM_2 2 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb SDM_10 10 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/ocf2_nocen_chr2_10kb SDM_20 20 1 out
########

#####BCF2
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb SDM_0 0 1 back ###filtering step off
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb SDM_1 1 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb SDM_2 2 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb SDM_10 10 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/bcf2_nocen_chr3_10kb SDM_20 20 1 back
########

#####mob1
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb SDM_0 0 1 back ###filtering step off
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb 0807_SDM_50 1 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb 0807_SDM_50 2 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb 0807_SDM_10 10 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/B_nocen_chr5_10kb 0807_SDM_5 20 1 back
########

#####mob2
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb SDM_0 0 1 back ###filtering step off
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb SDM_1 1 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb SDM_2 2 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb SDM_10 10 1 back
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/C_nocen_chr5_10kb SDM_20 20 1 back
########

#####sup1
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb SDM_0 0 1 out ###filtering step off
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb SDM_1 1 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb SDM_2 2 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb SDM_10 10 1 out
ruby SNP_distribution_method_variation.rb arabidopsis_datasets/No_centromere/sup1_nocen_chr4_10kb SDM_20 20 1 out
########

#####Analyse effect of ratio in model genome
ruby SNP_distribution_method_variation.rb Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left Ratio_0_1 0 1 back ##filtering step off
ruby SNP_distribution_method_variation.rb Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left Ratio_1_1 1 back
ruby SNP_distribution_method_variation.rb Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left Ratio_2_1 2 1 back
ruby SNP_distribution_method_variation.rb Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left Ratio_10_1 10 1 back
ruby SNP_distribution_method_variation.rb Small_genomes/arabidopsis_datasets/Analyse_effect_ratio/chr1_left Ratio_20_1 20 1 back
########


