#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'

require 'pp'
require 'benchmark'
require 'csv'

dataset = ARGV[0] 
# perm = ARGV[1]
# threshold = ARGV[2]

######Files
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
fasta_shuffle = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"


#Create list of fragment ids that contain  SNPs from vcf file 
snp_data, hm, ht, frag_pos = Stuff.snps_in_vcf(vcf_file)
# snp_data = ReformRatio.get_snp_data(vcf_file)




#Create hashes for fragments ids and SNP position

hm_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/ht_snps.txt")


# dic_pos_hm, dic_pos_ht =  Stuff.dic_id_pos(hm, ht, hm_list, ht_list)

# puts dic_pos_hm

# class Hash
#   def safe_invert
#     self.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
#   end
# end

# dic_pos_hm = dic_pos_hm.safe_invert

# dic_pos_ht = dic_pos_ht.safe_invert

dic_pos_hm = frag_pos[:hom]
dic_pos_ht = frag_pos[:het]

pp hm 

exit


##Create dictionaries with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm, dic_ht = Stuff.create_hash_snps(hm, ht)

##Create array with ordered fragments (fromf fasta_file) and from shuffled fragments (fasta_shuffle)
frags = ReformRatio.fasta_array(fasta_file)
frags_shuffled = ReformRatio.fasta_array(fasta_shuffle)

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok = ReformRatio.fasta_id_n_lengths(frags)
ids, lengths = ReformRatio.fasta_id_n_lengths(frags_shuffled)
genome_length = ReformRatio.genome_length(fasta_file)

##Define snps in hashes (fragment id as key and snp density as value). Create also lists 

##Assign the number of SNPs to each fragment in the shuffled list (hash)
##If a fragment does not have SNPs, the value assigned will be 0.
ok_hm, ok_ht, snps_hm, snps_ht = Stuff.define_snps(ids_ok, dic_hm, dic_ht)

#ratios
dic_ratios, ratios, ids_short = Stuff.important_ratios(snps_hm, snps_ht, ids_ok, threshold)

ids_short.flatten!


s_hm, s_ht, s_snps_hm, s_snps_ht = Stuff.define_snps(ids_short, dic_hm, dic_ht)

##if the ratio is not higher than a given threshold, eliminate the ids from the ids array
##and then eliminate those SNPs positions from the files, create new files.  
shuf_short_ids, hm_sh, ht_sh = Stuff.important_positions(ids_short, dic_pos_hm, dic_pos_ht, ids)

#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length
# dic_hm_norm, dic_ht_norm = Stuff.normalise_by_length(ids_, dic_hm, dic_ht, lengths)
dic_shuf_hm_norm, dic_shuf_ht_norm = Stuff.normalise_by_length(shuf_short_ids, dic_hm, dic_ht, lengths)

#Invert the hash so we can have the SNP density as a key.

dic_hm_inv = dic_shuf_hm_norm.safe_invert

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments 
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used 
#to reconstruct the right side of the distribution, and left, for the left side)

perm_hm = SDM.sorting(dic_hm_inv)

##Measuree time of SDM. Eventually add time needed for the remaining steps until we define the mutation
 Benchmark.bm do |b|
    b.report {10.times do ; perm_hm = SDM.sorting(dic_hm_inv);  end}
end

#Define SNPs in the recently ordered array of fragments.
dic_or_hm, dic_or_ht, snps_hm_or, snps_ht_or = Stuff.define_snps(perm_hm, dic_hm, dic_ht)

###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
#dic_ratios, ratios = Stuff.important_ratios(snps_hm, snps_ht, ids_ok)
dic_expected_ratios, expected_ratios = Stuff.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold)

#Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(perm_hm, frags_shuffled, ids)

#Create new fasta file with the ordered elements
File.open("arabidopsis_datasets/#{dataset}/frags_ordered#{threshold}.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered#{threshold}.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)

#Create arrays with the lists of SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)


###Calculate size of the group of fragments that have a high hm/ht ratio
contig_size = (genome_length/ids_ok.length).to_f
center = contig_size*(perm_hm.length)
puts "The length of the group of contigs that have a high hm/ht ratio is #{center.to_i} bp"


causal, candidate, percent = Mutation.define(hm_sh, ht_sh, hom_snps, het_snps, genome_length, ratios, expected_ratios)

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

distribution_plots = Mutation.distribution_plot(center, ratios, expected_ratios, dataset, perm)
