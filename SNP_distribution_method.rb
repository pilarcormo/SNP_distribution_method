#encoding: utf-8

require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/stuff'
require_relative 'lib/mutation'
require_relative 'lib/SDM'

require 'Bio'
require 'pp'
require 'csv'
require 'benchmark'
require 'PDist'

dataset = ARGV[0] # Name of dataset directory in 'small_genomes_SNPs/arabidopsis_datasets'
perm_files = ARGV[1]

######Files
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
fasta_shuffle = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"


#Create lists of SNPs
hm, ht = Stuff.snps_in_vcf(vcf_file)

##Create dictionaries with the id of the fragment as key and the NUMBER of SNP as value
dic_hm, dic_ht = Stuff.create_hash_snps(hm, ht)

##Open the fasta file with the randomly ordered fragments  and create an array with all the information
frags_shuffled = ReformRatio.fasta_array(fasta_shuffle)
frags = ReformRatio.fasta_array(fasta_file)

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok = ReformRatio.fasta_id_n_lengths(frags)
ids, lengths = ReformRatio.fasta_id_n_lengths(frags_shuffled)

ok_hm, ok_ht, snps_hm, snps_ht = Stuff.define_snps(ids_ok, dic_hm, dic_ht)

hm, ht, snps_hm_sh, snps_ht_sh = Stuff.define_snps(ids, dic_hm, dic_ht)

genome_length = ReformRatio.genome_length(fasta_file)

###Calculate ratios and delete those equal to or lower than 1 so only the important contigs remain.
x = 0
dic_ratios, ratios = {}, []
snps_hm.length.times do
	ratio = (snps_hm[x]+1)/(snps_ht[x]+1)
	dic_ratios.store(ids_ok[x], ratio.to_f) 
	x = x + 1
end

dic_ratios.delete_if { |id, ratio|  ratio <= 1  }

ratios << dic_ratios.values
ratios.flatten!

###Use divisions
# ratios_div = ratios.each_slice(div).to_a

# sum = []
# ratios_div.each do |array|
# 	sum << array.inject(0) {|sum, i|  sum + i }
# end


# ratios.each_slice(div) do |x,y,z|
# 	ratios_div << (x+y+z).to_f
# end 

# pp ratios_div

# ratios.each do |array|
# 	array = array[0].to_f + array[1].to_f + array[2].to_f
# 	pp array
# end

# pp ratios




##Assign the number of SNPs to each fragment in the shuffled list. 
##If a fragment does not have SNPs, the value assigned will be 0.

#

#Define SNPs per fragment in the shuffled fasta array and then normalise the value of SNP density per fragment length

dic_shuf_hm_norm, dic_shuf_ht_norm = Stuff.normalise_by_length(ids, dic_hm, dic_ht, lengths)
dic_hm_norm, dic_ht_norm = Stuff.normalise_by_length(ids_ok, dic_hm, dic_ht, lengths)

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

dic_or_hm, dic_or_ht, snps_hm_or, snps_ht_or = Stuff.define_snps(perm_hm, dic_hm, dic_ht)


x = 0
dic_expected_ratios = {}
expected_ratios = []
snps_hm_or.length.times do
	expected_ratio = (snps_hm_or[x]+1)/(snps_ht_or[x]+1)
	dic_expected_ratios.store(perm_hm[x], expected_ratio.to_f)
	x = x + 1
end

#Remove ratios when they are lower than or equal to 1. 
deleted_ids = []
#revise this, maybe 1 is to restrictive, just to try for now 

dic_expected_ratios.each do |id, ratio|
	if ratio <= 1
		deleted_ids << id
	end
end 
# expected_ratios << dic_expected_ratios.values

deleted_ids.each { |element|
  if ids_ok.include?(element)
    ids_ok.delete(element)
  end
}

ok_hm_del, ok_ht_del, snps_hm_del, snps_ht_del = Stuff.define_snps(ids_ok, dic_hm, dic_ht)

dic_expected_ratios.delete_if { |id, ratio|  ratio <= 1 } 
short_sh = []

short_sh = dic_expected_ratios.keys

ratios.shift
ratios.flatten!
expected_ratios.flatten!

hm_del, ht_del, snps_hm_del2, snps_ht_del2 = Stuff.define_snps(short_sh, dic_hm, dic_ht)

##Take IDs, lenght and sequence from the shuffled fasta file and add them to the permutation array 

fasta_perm = Stuff.create_perm_fasta(short_sh, frags_shuffled, ids)
# fasta_ok = Stuff.create_perm_fasta(short_ok, frags, ids_ok)


#Create new fasta file with the ordered elements
File.open("arabidopsis_datasets/#{dataset}/frags_ordered.fasta", "w+") do |f|
  fasta_perm.each { |element| f.puts(element) }
end

fasta_ordered = "arabidopsis_datasets/#{dataset}/frags_ordered.fasta"
frags_ordered = ReformRatio.fasta_array(fasta_ordered)

#Create arrays with the lists of SNP positions in the new ordered file.

snp_data = ReformRatio.get_snp_data(vcf_file)

het_snps, hom_snps = ReformRatio.perm_pos(frags_ordered, snp_data)

###Calculated size of the group of fragments that have a high hm to ht ratio

contig_size = (genome_length/perm_hm.length).to_f
center = contig_size*(short_sh.length)

###Create arrays of correct SNP positions
hm_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/ht_snps.txt")

hm_list_2, hm_list_3 = [], []
hm_list_2 << hm_list
hm_list_3 << hm_list

hm_list_2.flatten!
hm_list_3.flatten!

pp ok_hm_del

hm_ok = Stuff.positions_by_fragment(ok_hm_del, hm_list)

hm_p = Stuff.positions_by_fragment(hm_del, hm_list_2)




####Eventually do a new method in Stuff with this
##Create a hash with the fragments ids as keys and the SNP positions per fragment as value 

positions_hm = []
short_sh.each do |frag|
	if hm_p.has_key?(frag)
		positions_hm << hm_p[frag]
	end 
end
positions_hm.flatten!

# mutation = Mutation.define(hm_list_3, ht_list, positions_hm, het_snps, genome_length, ratios, expected_ratios)

# puts positions_hm

Dir.mkdir("arabidopsis_datasets/#{dataset}/#{perm_files}")
Dir.chdir("arabidopsis_datasets/#{dataset}/#{perm_files}") do
	WriteIt::write_txt("perm_hm", positions_hm) # save the SNP distributions for the best permutation in the generation
	WriteIt::write_txt("perm_ht", het_snps)
end

##Measuree time of SDM. Eventually add time needed for the remaining steps until we define the mutation

 Benchmark.bm do |b|
    b.report {10.times do ; perm_hm = SDM.sorting(dic_hm_inv);  end}
end

