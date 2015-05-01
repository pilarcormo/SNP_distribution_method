#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'
require_relative 'lib/plot'

require 'pp'
require 'benchmark'
require 'csv'


if ARGV.empty?
	puts "Please specify a (1) dataset, a (2) name for the output folder, a (3) threshold to discard the contigs which a ratio below it and a (4) factor to calculate the ratio (1, 0.1, 0.01...) "
else 
	dataset = ARGV[0] 
	perm = ARGV[1]
	threshold = ARGV[2].to_i
	adjust = ARGV[3]
	puts "Looking for SNPs in #{dataset}"
	puts "Output will be in #{dataset}/#{perm}"
	puts "A factor of #{adjust} will be used to calculate the ratio"
	if threshold == 1
		puts "Filtering step on"
	elsif threshold == 0
		puts "Filtering step off"
	else 
		puts "Not valid filtering value, plese specify 0 to skip filtering and 1 to allow it"
		exit
	end 
	#puts "All the contigs with a ratio lower than #{threshold} will not be considered in the analysis"
end 

###############################################################
######Files
loc = "arabidopsis_datasets/#{dataset}"
vcf_file = "#{loc}/snps.vcf"
fasta_file = "#{loc}/frags.fasta"
fasta_shuffle = "#{loc}/frags_shuffled.fasta"

hm_list = WriteIt.file_to_ints_array("#{loc}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("#{loc}/ht_snps.txt")
##############################################################
##############################################################

#Create arrays of ids of those fragments that contain  SNPs in the vcf file and hashes with the id and position eh
snp_data, hm, ht, frag_pos = Stuff.snps_in_vcf(vcf_file)

frag_pos_hm = frag_pos[:hom]
frag_pos_ht = frag_pos[:het]

#Create hashes for fragments ids and SNP position for the correct genome - reference 

dic_pos_hm =  Stuff.dic_id_pos(hm, hm_list)
dic_pos_ht =  Stuff.dic_id_pos(ht, ht_list)
 ########

##Create dictionaries with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm = Stuff.create_hash_number(hm)
dic_ht = Stuff.create_hash_number(ht)

##Create array with ordered fragments (fromf fasta_file) and from shuffled fragments (fasta_shuffle)
frags = ReformRatio.fasta_array(fasta_file)
frags_shuffled = ReformRatio.fasta_array(fasta_shuffle)

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok, id_len_ok = ReformRatio.fasta_id_n_lengths(frags)
ids, lengths, id_len = ReformRatio.fasta_id_n_lengths(frags_shuffled)

genome_length = ReformRatio.genome_length(fasta_file)

##Define snps in hashes (fragment id as key and snp density as value). Create also lists 

##Assign the number of SNPs to each fragment in the shuffled list (hash)
##If a fragment does not have SNPs, the value assigned will be 0.
ok_hm, snps_hm = Stuff.define_snps(ids_ok, dic_hm)
ok_ht, snps_ht = Stuff.define_snps(ids_ok, dic_ht)

#ratios

dic_ratios, ratios, ids_short = Stuff.important_ratios(snps_hm, snps_ht, ids_ok, threshold, adjust)


s_hm, s_snps_hm = Stuff.define_snps(ids_short, dic_hm)
s_ht, s_snps_ht = Stuff.define_snps(ids_short, dic_ht)

##if the ratio is not higher than a given threshold, eliminate the ids from the ids array
##and then eliminate those SNPs positions from the files, create new files.  
shuf_short_ids = Stuff.important_ids(ids_short, ids)
hm_sh = Stuff.important_pos(ids_short, dic_pos_hm)
ht_sh = Stuff.important_pos(ids_short, dic_pos_ht)


shuf_hm, shu_snps_hm = Stuff.define_snps(shuf_short_ids, dic_hm)


#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length

dic_shuf_hm_norm = Stuff.normalise_by_length(lengths, shuf_hm)

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments 
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used 
#to reconstruct the right side of the distribution, and left, for the left side)
puts "\n"
perm_hm, mut = SDM.sorting(dic_shuf_hm_norm)

perm_hm

# snps_around_mutation = Stuff.important_pos(mut, dic_pos_hm)

# pp snps_around_mutation

##Measuree time of SDM. Eventually add time needed for the remaining steps until we define the mutation
puts "Time spent sorting the contigs:"
Benchmark.bm do |b|
    b.report {10.times do ; perm_hm = SDM.sorting(dic_shuf_hm_norm);  end}
end
puts "done"

#Define SNPs in the recently ordered array of fragments.
dic_or_hm, snps_hm_or = Stuff.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Stuff.define_snps(perm_hm, dic_ht)

###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
#dic_ratios, ratios = Stuff.important_ratios(snps_hm, snps_ht, ids_ok)
dic_expected_ratios, expected_ratios, ids_short = Stuff.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold, adjust)

csv = "#{loc}/ratio_positions#{threshold}_#{adjust}.csv"
csv_perm = "#{loc}/ratio_positions_perm_#{threshold}_#{adjust}.csv"
pos_ratios = Stuff.csv_pos_ratio(csv, dic_pos_hm, dic_ratios)
pos_ratios_perm = Stuff.csv_pos_ratio(csv_perm, dic_pos_hm, dic_expected_ratios)
#Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(perm_hm, frags_shuffled, ids)

#Create new fasta file with the ordered elements
File.open("#{loc}/frags_ordered_thres#{threshold}.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered#{threshold}.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)
ids_or, lengths_or, id_len_or = ReformRatio.fasta_id_n_lengths(frags_ordered)


dic_global_pos, hm_global_positions_perm = Stuff.define_global_pos(perm_hm, frag_pos_hm, id_len_or) 
dic_global_pos, ht_global_positions_perm = Stuff.define_global_pos(perm_hm, frag_pos_ht, id_len_ok) 


WriteIt::write_txt("arabidopsis_datasets/#{dataset}/global_perm_hm", hm_global_positions_perm)
expected = WriteIt.file_to_ints_array("Reads/Aw_sup1-2/Variant_calling/sup1_2_1/hm_nocen.txt")
Plot::exp_vs_hyp_densities(expected, hm_global_positions_perm, loc)
Plot::densities(hm_global_positions_perm, ht_global_positions_perm, pos_ratios_perm, genome_length, loc)
exit 


#Create arrays with the lists of SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)
puts "\n"
###Calculate size of the group of fragments that have a high hm/ht ratio
contig_size = (genome_length/ids_ok.length).to_f
center = contig_size*(perm_hm.length)
puts "The length of the group of contigs that have a high hm/ht ratio is #{center.to_i} bp"
puts "..."
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)


"Defining the mutation..."

causal, candidate, percent = Mutation.define(hm_sh, ht_sh, hom_snps, het_snps, genome_length, ratios, expected_ratios, pos_ratios, pos_ratios_perm)
puts "\n"
puts "Writing the output..."
Dir.mkdir("arabidopsis_datasets/#{dataset}/#{perm}")
Dir.chdir("arabidopsis_datasets/#{dataset}/#{perm}") do
	WriteIt::write_txt("perm_hm", hom_snps) # save the SNP distributions for the best permutation in the generation
	WriteIt::write_txt("perm_ht", het_snps)
  WriteIt::write_txt("hm_snps_short", hm_sh) # save the SNP distributions for the best permutation in the generation
  WriteIt::write_txt("ht_snps_short", ht_sh)
	File.open("mutation.txt", "w+") do |f|
		f.puts "The length of the group of contigs that form the peak of the distribution is #{center.to_i} bp"
		f.puts "Location of causal mutation in correctly ordered genome: #{causal}"
		f.puts "Candidate SNP position in permutation: #{candidate}"
		f.puts "Shift #{percent} %"
	end
end

distribution_plots = Mutation.distribution_plot(center, ratios, expected_ratios, dataset, perm, pos_ratios_perm)


