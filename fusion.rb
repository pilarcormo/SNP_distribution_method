
#encoding: utf-8
require_relative 'lib/write_it'
require 'pp'

location = ARGV[0]

chr1 = 34964571
chr2 = 22037565
chr3 = 25499034
chr4 = 20862711

hm_list1 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome1/interesting_1/hm_nocen.txt") 
ht_list1 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome1/interesting_1/ht_nocen.txt")

hm_list2 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome2/interesting_2/hm_nocen.txt") 
ht_list2 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome2/interesting_2/ht_nocen.txt")

hm_list3 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome3/interesting_3/hm_nocen.txt") 
ht_list3 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome3/interesting_3/ht_nocen.txt")

hm_list4 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome4/interesting_4/hm_nocen.txt") 
ht_list4 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome4/interesting_4/ht_nocen.txt")

hm_list5 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome5/interesting_5/hm_nocen.txt") 
ht_list5 = WriteIt.file_to_ints_array("Reads/#{location}/C_chromosome5/interesting_5/ht_nocen.txt")


class Genome 
	def self.create_whole(snps, chr_length)
		whole_size = []
		snps.each do |position|
			new_pos = position.to_i + chr_length.to_i
			whole_size << new_pos
		end 
		return whole_size
	end 
end 
chr12 = chr1 + chr2
chr123 = chr12 + chr3
chr1234 = chr123 + chr4

chrtotal = chr1234 + chr5
puts chrtotal
exit 


hm_nocen2 = Genome.create_whole(hm_list2, chr1)
ht_nocen2 = Genome.create_whole(hm_list2, chr1)
hm_nocen3 = Genome.create_whole(hm_list3, chr12)
ht_nocen3 = Genome.create_whole(hm_list3, chr12)
hm_nocen4 = Genome.create_whole(hm_list4, chr123)
ht_nocen4 = Genome.create_whole(hm_list4, chr123)
hm_nocen5 = Genome.create_whole(hm_list5, chr1234)
ht_nocen5 = Genome.create_whole(hm_list5, chr1234)

hm_nocen_whole, ht_nocen_whole = [], []

hm_nocen_whole << hm_list1
hm_nocen_whole << hm_nocen2
hm_nocen_whole << hm_nocen3
hm_nocen_whole << hm_nocen4
hm_nocen_whole << hm_nocen5
ht_nocen_whole << ht_list1
ht_nocen_whole << ht_nocen2
ht_nocen_whole << ht_nocen3
ht_nocen_whole << ht_nocen4
ht_nocen_whole <<  ht_nocen5

hm_nocen_whole.flatten!
ht_nocen_whole.flatten!

final = "Reads/#{location}/C_wholegenome"
Dir.mkdir("#{final}")

File.open("#{final}/hm_nocen_whole.txt", "w+") do |f|
  hm_nocen_whole.each { |element| f.puts(element) }
end

File.open("#{final}/ht_nocen_whole.txt", "w+") do |f|
  ht_nocen_whole.each { |element| f.puts(element) }
end