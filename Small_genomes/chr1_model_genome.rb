#encoding: utf-8
require_relative 'lib/model_genome'
require_relative '../lib/write_it'
require_relative  'lib/reform_ratio'
require 'pp'

name = ARGV[0]
contig_size = ARGV[1].to_i

fasta_file = "TAIR10_chr1.fasta"

genome_length = ReformRatio::genome_length(fasta_file)

genome_length.to_i

snp = genome_length/3000
snp_2 = snp*2
size_2 = genome_length/2
# make the directory to put data files into
Dir.mkdir("arabidopsis_datasets/#{name}")

hm_r = "hm <- rnorm(#{snp}, #{size_2}, #{snp_2})" # Causative SNP at/near 10000
ht_r = "ht <- runif(#{snp}, 1, #{genome_length})"   # Genome length of 10000
hm, ht = ModelGenome::get_snps(hm_r, ht_r)

arabidopsis_c1 = ModelGenome::fasta_to_char_array(fasta_file)

puts "There are #{hm.length} homozygous SNPs"
puts "There are #{ht.length} heterozygous SNPs"

snp_pos = [hm, ht].flatten

puts "...and generating the fragments"

frags = ModelGenome::get_frags(arabidopsis_c1, contig_size)

puts "Small genome length: #{genome_length} bases"
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle
vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_data("arabidopsis_datasets/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("arabidopsis_datasets/#{name}/snps.vcf", vcf)
WriteIt::write_data("arabidopsis_datasets/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)
WriteIt::write_txt("arabidopsis_datasets/#{name}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
WriteIt::write_txt("arabidopsis_datasets/#{name}/hm_snps", hm)
WriteIt::write_txt("arabidopsis_datasets/#{name}/ht_snps", ht)