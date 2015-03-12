#encoding: utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'
require_relative 'lib/reform_ratio'
require 'pp'

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "SNP_distribution_method/arabidopsis_datasets/#{ARGV[0]}"))

# Create the lists of homozygous and heterozygous SNPs
fasta_file = "TAIR10_Chr1.fasta"
genome_length = ReformRatio::genome_length(fasta_file)

# hm_r = "hm <- rnorm(10000, 15213836, 1000000)" # Causative SNP at/near 1000
# ht_r = "ht <- runif(10000, 1, #{genome_length})" # Genome length is same as arabidopisis chromosome 4
# hm, ht = ModelGenome::get_snps(hm_r, ht_r)
# snp_pos = [hm, ht].flatten

hm = WriteIt::file_to_array("chromosome1/hm.txt")
ht = WriteIt::file_to_array("chromosome1/ht.txt")

hm_r = "hm <- rnorm(5000, 10, 1000000)" 
ht_r = "ht <- runif(10000, 1, #{genome_length})"
hm2, ht2 = ModelGenome::get_snps(hm_r, ht_r)

pp hm
pp ht
pp hm2
pp ht2

# puts "There are #{hm.length} homozygous SNPs"
# puts "There are #{ht.length} heterozygous SNPs"
# puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?(15213836)}"

# arabidopsis_c1 = ModelGenome::fasta_to_char_array(fasta_file)
# puts "Fragmenting arabidopsis chromosome 1..."
# contig_size = 10000 # 10-20kb
# frags = ModelGenome::get_frags(arabidopsis_c1, contig_size)
# puts "Done!"
# puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"


# # Get the positions of the SNPs on fragments
# pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

# fastaformat_array = ModelGenome::fasta_array(frags)
# fastaformat_array_shuf = fastaformat_array.shuffle # shuffle it to show that the order doesn't need to be conserved when working out density later on

# vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

# WriteIt::write_data("chromosome1/#{ARGV[0]}/frags.fasta", fastaformat_array)
# WriteIt::write_data("chromosome1/#{ARGV[0]}/snps.vcf", vcf)
# WriteIt::write_data("chromosome1/#{ARGV[0]}/frags_shuffled.fasta", fastaformat_array_shuf)
# WriteIt::write_txt("chromosome1/#{ARGV[0]}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
WriteIt::write_txt("chromosome1/#{ARGV[0]}/hm_snps", hm)
WriteIt::write_txt("chromosome1/#{ARGV[0]}/ht_snps", ht)