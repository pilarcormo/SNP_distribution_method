#encoding: utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'
require_relative 'lib/reform_ratio'

name = ARGV[0]
contig_size = ARGV[1].to_i
file = ARGV[2]
chr = ARGV[3]

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "SNP_distribution_method/arabidopsis_datasets/#{name}"))

# Create the lists of homozygous and heterozygous SNPs
fasta_file = "TAIR10_chr#{chr}.fasta"
genome_length = ReformRatio::genome_length(fasta_file)

#####Real hm and ht SNPs 
hm = WriteIt::file_to_ints_array("#{file}/hm_nocen.txt")
ht = WriteIt::file_to_ints_array("#{file}/ht_nocen.txt")


snp_pos = [hm, ht].flatten

puts "There are #{hm.length} homozygous SNPs"
puts "There are #{ht.length} heterozygous SNPs"


arabidopsis_c = ModelGenome::fasta_to_char_array(fasta_file)
puts "Fragmenting arabidopsis chromosome #{chr}..."
frags = ModelGenome::get_frags(arabidopsis_c1, contig_size)
puts "Done!"
puts "Arabidopsis chr length: #{arabidopsis_c.length} kb"
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle
vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_txt("arabidopsis_datasets/#{name}/hm_snps", hm)
WriteIt::write_txt("arabidopsis_datasets/#{name}/ht_snps", ht)
WriteIt::write_data("arabidopsis_datasets/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("arabidopsis_datasets/#{name}/snps.vcf", vcf)
WriteIt::write_data("arabidopsis_datasets/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)

