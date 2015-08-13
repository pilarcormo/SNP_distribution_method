require 'bio'
require 'bio-samtools'	
require 'pp'
# class Contigs


# 	def self.fasta_array(fasta_file)

# 		return fasta
# 	end

# 	def self.fasta_id_n_lengths(fasta)
# 		ids, lengths = [], []
# 		id_len = {}
# 		fasta.each do |i|
# 			ids << i.entry_id
# 			lengths << i.length
# 			id_len.store(i.entry_id, i.length)
# 		end
# 		return ids, lengths, id_len
# 	end
# end 


fasta_file = "contigs.fasta"

fasta = [] # we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
Bio::FastaFormat.open(fasta_file).each do |i| # get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
	fasta << i
end

ids, lengths = [], []
id_len = {}
fasta.each do |i|
	ids << i.entry_id
end
ids.each do |i|
	a =  i.split("_")
	lengths << a[3]
end 

pp lengths