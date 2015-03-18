#encoding: utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'
require_relative 'lib/reform_ratio'

name = ARGV[0]
contig_size = ARGV[1].to_i
chr = ARGV[2]

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "Pilar/SDM/genomes/#{name}"))

# Create the lists of homozygous and heterozygous SNPs
fasta_file = "TAIR10_chr#{chr}.fasta"
genome_length = ReformRatio::genome_length(fasta_file)

#####Real hm and ht SNPs 
hm = WriteIt::file_to_ints_array("BCF2_chromosome#{chr}/hm.txt")
ht = WriteIt::file_to_ints_array("BCF2_chromosome#{chr}/ht.txt")

# hm_r = "hm <- rnorm(5000, #{snp}, 1000000)" 
# ht_r = "ht <- runif(10000, 1, #{genome_length})"
# hm, ht = ModelGenome::get_snps(hm_r, ht_r)


snp_pos = [hm, ht].flatten

# puts "There are #{hm.length} homozygous SNPs"
# puts "There are #{ht.length} heterozygous SNPs"
# puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?()}"

arabidopsis_c = ModelGenome::fasta_to_char_array(fasta_file)
puts "Fragmenting arabidopsis chromosome 1..."
# 10-20kb
frags = ModelGenome::get_frags(arabidopsis_c1, contig_size)
puts "Done!"
puts "Arabidopsis chr length: #{arabidopsis_c.length}/ kb"
puts "Fragmented seq   length: #{frags.join.length} = close enough? You decide."
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle
vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_txt("genomes/#{name}/hm_snps", hm)
WriteIt::write_txt("genomes/#{name}/ht_snps", ht)
WriteIt::write_data("genomes/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("genomes/#{name}/snps.vcf", vcf)
WriteIt::write_data("genomes/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)
# WriteIt::write_txt("genomes/#{name}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
