#encoding: utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'
require 'pp'


name = ARGV[0]
size = 2000

# snp = (size/1000)*2

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "SNP_distribution_method/arabidopsis_datasets/#{name}"))

# Create the lists of homozygous and heterozygous SNPs
# hm_r = 'hm <- rnorm(#{snp}, #{size}/2, #{snp}*2)' # Causative SNP at/near 10000
# ht_r = 'ht <- runif(#{snp}, 1, #{size})'   # Genome length of 10000
# hm, ht = ModelGenome::get_snps(hm_r, ht_r)


# puts "There are #{hm.length} homozygous SNPs"
# puts "There are #{ht.length} heterozygous SNPs"
# puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?(7500000)}"
hm, ht = [], []

hm = WriteIt::file_to_ints_array("chromosom1/hm.txt")
ht = WriteIt::file_to_ints_array("chromosome/ht.txt")
# 

# hm = WriteIt::file_to_array("hm.txt")
# ht = WriteIt::file_to_array("ht.txt")



# hm = [1000, 1001, 1003, 1004, 1010, 1020]
# ht = [25, 125, 230, 856, 965, 1100, 1229, 1225, 1436]

# hm.flatten!
# ht.flatten!

snp_pos = [hm, ht].flatten


arabidopsis_c4 = ModelGenome::fasta_to_char_array("TAIR10_chr4.fasta")
puts "Creating the genome..."
small_genome = arabidopsis_c4[-size..-1] # Genome length of 100 kb
contig_size = 10 # 100-200 bp
puts "...and generating the fragments"
frags = ModelGenome::get_frags(small_genome, contig_size)

puts "Small genome     length: #{small_genome.length} bases"
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"


# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)


fastaformat_array = ModelGenome::fasta_array(frags)
vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_data("arabidopsis_datasets/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("arabidopsis_datasets/#{name}/snps.vcf", vcf)
# WriteIt::write_data("arabidopsis_datasets/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)
# WriteIt::write_txt("arabidopsis_datasets/#{name}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
WriteIt::write_txt("arabidopsis_datasets/#{name}/hm_snps", hm)
WriteIt::write_txt("arabidopsis_datasets/#{name}/ht_snps", ht)