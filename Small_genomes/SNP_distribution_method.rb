#encoding: utf-8
require_relative '../lib/reform_ratio'
require_relative '../lib/write_it'
require_relative '../lib/reform_ratio'
require_relative '../lib/stuff'
require_relative 'lib/mutation'
require_relative '../lib/SDM'
require_relative '../lib/ratio_filtering'
require_relative '../lib/output'
require_relative 'lib/snp_dist'

require 'pp'
require 'benchmark'
require 'csv'


if ARGV.empty?
  puts "Please specify a (1) dataset, a (2) name for the output folder, a (3) threshold to discard the contigs which a ratio below it"
else 
  dataset = ARGV[0] 
  output_folder = ARGV[1]
  threshold = ARGV[2].to_i #degree of filtering:  100, 50, 10, 5
  adjust = ARGV[3]
  cross = ARGV[4]
  puts "Looking for SNPs in #{dataset}"
  puts "Output will be in #{dataset}/#{output_folder}"
  if threshold > 0
    puts "Filtering step on: #{threshold}% selected"
  elsif threshold == 0
    puts "Filtering step off. "
  else 
    puts "Not valid filtering value, plese specify 0 to skip filtering or a positive integer to allow it"
    exit
  end 
  puts "A factor of #{adjust} will be used to calculate the ratio"
  puts "#{cross}-cross selected"
end 

######Files
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
fasta_shuffle = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"
hm_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/hm_snps.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/ht_snps.txt")
##############################################################
##############################################################


####[1] Open VCF file 
snp_data, hm, ht, frag_pos = Stuff.snps_in_vcf(vcf_file)

##Hash with the fragments ids and the homozygous and heterozygous SNP positions
frag_pos_hm = frag_pos[:hom] 
frag_pos_ht = frag_pos[:het]

##Hashes with fragments ids and SNP positions for the correctly ordered genome
dic_pos_hm =  Stuff.dic_id_pos(hm, hm_list)
dic_pos_ht =  Stuff.dic_id_pos(ht, ht_list)


##Create hashes with the id of the fragment as the key and the NUMBER of SNPs as value
dic_hm = Stuff.create_hash_number(hm)
dic_ht = Stuff.create_hash_number(ht)


####[2] Open FASTA files containing the ordered and unordered contigs
##Create array with ordered fragments (from fasta_file) and from shuffled fragments (fasta_shuffle)
frags = ReformRatio.fasta_array(fasta_file)
frags_shuffled = ReformRatio.fasta_array(fasta_shuffle)

##From the previous array take ids and lengths and put them in 2 separate new arrays
ids_ok, lengths_ok, id_len_ok = ReformRatio.fasta_id_n_lengths(frags)
ids, lengths, id_len = ReformRatio.fasta_id_n_lengths(frags_shuffled)

genome_length = ReformRatio.genome_length(fasta_file)

##Assign the number of SNPs to each fragment in the shuffled list (hash)
##If a fragment does not have SNPs, the value assigned will be 0.
ok_hm, snps_hm = Stuff.define_snps(ids_ok, dic_hm)
ok_ht, snps_ht = Stuff.define_snps(ids_ok, dic_ht)
shuf_hm, shuf_snps_hm = Stuff.define_snps(ids, dic_hm)
shuf_ht, shuf_snps_ht = Stuff.define_snps(ids, dic_ht)


####[3] Pre-filtering. Calculate Hom/het ratios for shuffled and ordered set of contigs. If a threshold > 0 was provided, this value is the percentage of the maximum ratio below which a contig will be discarded

dic_ratios, ratios, ids_short, dic_ratios_inv  = Ratio_filtering.important_ratios(snps_hm, snps_ht, ids_ok, threshold, adjust) 
dic_ratios_shuf, ratios_shuf, ids_short_shuf, dic_ratios_inv_shuf = Ratio_filtering.important_ratios(shuf_snps_hm, shuf_snps_ht, ids, threshold, adjust) 

##Redefine the arrays of SNPs after discarding the contigs that fell below the threshold provided. We refered to them as the "important contigs" and the SNPs on those are the "important positions"

s_hm, s_snps_hm = Stuff.define_snps(ids_short, dic_hm)

shuf_short_ids = Ratio_filtering.important_ids(ids_short, ids)
hm_sh = Ratio_filtering.important_pos(ids_short, dic_pos_hm)
ht_sh = Ratio_filtering.important_pos(ids_short, dic_pos_ht)

shuf_hm, shu_snps_hm = Stuff.define_snps(shuf_short_ids, dic_hm)

##Calculate how many contigs were discarded

contigs_discarded = ids.length - ids_short.length
puts "#{contigs_discarded} contigs out of #{ids.length} discarded"

####[4] SDM 
## Calculate scores (number of homozygous SNPs in each contig divided by fragment length)

dic_shuf_hm_norm = SDM.normalise_by_length(lengths, shuf_hm)

##Iteration: look for the minimum value in the array of values, that will be 0 (fragments without SNPs) and put the fragments 
#with this value in a list. Then, the list is cut by half and each half is added to a new array (right, that will be used 
#to reconstruct the right side of the distribution, and left, for the left side)

perm_hm, perm_ratio, mut, hyp_positions = SDM.calling_SDM(dic_shuf_hm_norm, dic_ratios_inv_shuf, cross, dic_pos_hm)

puts "Hypothetical positions carrying the causal mutation #{hyp_positions}"


#Define SNPs in the r ordered array of fragments.
dic_or_hm, snps_hm_or = Stuff.define_snps(perm_hm, dic_hm)
dic_or_ht, snps_ht_or = Stuff.define_snps(perm_hm, dic_ht)

#Calculate ratios in the contig permutation obtained from SDM
dic_expected_ratios, expected_ratios, exp_ids_short, exp_inv_ratios = Ratio_filtering.important_ratios(snps_hm_or, snps_ht_or, perm_hm, threshold, adjust)

####[5] Outputs

#Create FASTA file for the contig permutation obtained from SDM
fasta_perm = Output.create_perm_fasta(perm_hm, frags_shuffled, ids)

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

causal, candidate, percent = Mutation.define(hm_sh, ht_sh, hom_snps, het_snps, center, ratios, expected_ratios)

Dir.mkdir("arabidopsis_datasets/#{dataset}/#{output_folder}")
Dir.chdir("arabidopsis_datasets/#{dataset}/#{output_folder}") do
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


if File.exist?('1-15Mb.csv') == true 
  CSV.open("1-15Mb.csv", "a+") do |csv|
    csv << [genome_length, ids_ok.length, percent]
  end 
else 
  CSV.open("1-15Mb.csv", "wb") do |csv|
    csv << ["Genome_size","Contigs","Deviation"]
  end 
  CSV.open("1-15Mb.csv", "a+") do |csv|
    csv << [genome_length, ids_ok.length, percent]
  end 
end 

####[6] Plots
distribution_plots = Mutation.distribution_plot(center, ratios, expected_ratios, dataset, output_folder)