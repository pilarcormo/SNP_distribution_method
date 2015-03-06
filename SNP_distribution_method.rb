#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'
require_relative 'lib/snp_dist'

require 'pp'
require 'benchmark'

dataset = ARGV[0] # Name of dataset directory in 'small_genomes_SNPs/arabidopsis_datasets'
perm = ARGV[1]

######Files
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
fasta_shuffle = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"


#Create lists of SNPs
hm, ht = Stuff.snps_in_vcf(vcf_file)
snp_data = ReformRatio.get_snp_data(vcf_file)

##Create dictionaries with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm, dic_ht = Stuff.create_hash_snps(hm, ht)

##Open the fasta file with the randomly ordered fragments  and create an array with all the information
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


#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length
dic_hm_norm, dic_ht_norm = Stuff.normalise_by_length(ids_ok, dic_hm, dic_ht, lengths)
dic_shuf_hm_norm, dic_shuf_ht_norm = Stuff.normalise_by_length(ids, dic_hm, dic_ht, lengths)
pp dic_hm_norm

#Invert the hash so we can have the SNP density as a key.

class Hash
  def safe_invert
    self.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
  end
end

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

#######NEEDED????? CHECK THIS:
dic_or_hm, dic_or_ht, snps_hm_or, snps_ht_or = Stuff.define_snps(perm_hm, dic_hm, dic_ht)

###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
dic_ratios, ratios = Stuff.important_ratios(snps_hm, snps_ht, ids_ok)
dic_expected_ratios, expected_ratios = Stuff.important_ratios(snps_hm_or, snps_ht_or, perm_hm)



##Create a shorter version of the ordered array of fragments with only the fragments that have a high hm/ht ratio 


short_ids = []
dic_expected_ratios.each do |id, ratio|
	short_ids << id
end 

ids_ok_short = []
ids_ok.each { |element|
	if short_ids.include?(element)
		ids_ok_short << element
	else 
    	ids_ok.delete(element)
  	end
}

ok_hm_del, ok_ht_del, snps_hm_del, snps_ht_del = Stuff.define_snps(ids_ok, dic_hm, dic_ht)

short_or = []
short_or = dic_expected_ratios.keys

hm_del, ht_del, snps_hm_del2, snps_ht_del2 = Stuff.define_snps(short_or, dic_hm, dic_ht)


##Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(perm_hm, frags_shuffled, ids)
fasta_perm_short = Stuff.create_perm_fasta(short_or, frags_shuffled, ids)

#Create new fasta file with the ordered elements
File.open("arabidopsis_datasets/#{dataset}/frags_ordered.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end
File.open("arabidopsis_datasets/#{dataset}/frags_ordered_short.fasta", "w+") do |f|
  fasta_perm_short.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)

#Create arrays with the lists of SNP positions in the new ordered file.
het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)

###Calculate size of the group of fragments that have a high hm/ht ratio
contig_size = (genome_length/perm_hm.length).to_f
center = contig_size*(short_or.length)
puts "The length of the group of contigs that have a high hm/ht ratio is #{center.to_i} bp"

###Create arrays of correct SNP positions
hm_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/ht_snps.txt")

hm_list_2, hm_list_3 = [], []
hm_list_2 << hm_list
hm_list_2.flatten!
hm_list_3 << hm_list
hm_list_3.flatten!


dic_positions_ok = Stuff.positions_by_fragment(ok_hm_del, hm_list)

dic_positions_or = Stuff.positions_by_fragment(hm_del, hm_list_2)


####Eventually do a new method in Stuff with this
##Create a hash with the fragments ids as keys and the SNP positions per fragment as value 

positions_hm = []
short_or.each do |frag|
	if dic_positions_or.has_key?(frag)
		positions_hm << dic_positions_or[frag]
	end 
end
positions_hm.flatten!

pp short_ids.length
pp short_or.length


causal, candidate, percent = Mutation.define(hm_list_3, ht_list, positions_hm, het_snps, genome_length, ratios, expected_ratios)


Dir.mkdir("arabidopsis_datasets/#{dataset}/#{perm}")
Dir.chdir("arabidopsis_datasets/#{dataset}/#{perm}") do
	WriteIt::write_txt("perm_hm", positions_hm) # save the SNP distributions for the best permutation in the generation
	WriteIt::write_txt("perm_ht", het_snps)
	File.open("mutation.txt", "w+") do |f|
		f.puts "The length of the group of contigs that form the peak of the distribution is #{center.to_i} bp"
		f.puts "Location of causal mutation in correctly ordered genome: #{causal}"
		f.puts "Candidate SNP position in permutation: #{candidate}"
		f.puts "Shift #{percent} %"
	end
end


distribution_plots = Mutation.distribution_plot(center, ratios, expected_ratios, dataset, perm)



