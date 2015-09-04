
#encoding: utf-8
require_relative 'lib/write_it'
require 'pp'

chromosome = ARGV[0] 
location = ARGV[1]

hm_list = WriteIt.file_to_ints_array("Reads/#{location}/hm.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("Reads/#{location}/ht.txt")

centromere = {"chr1" => [15086545-3950000, 15086545+3950000],"chr2" => [3608429- 1500000, 3608429+ 1500000], "chr3" => [14209452- 1500000, 14209452+1500000], "chr4" => [3956521- 1400000, 3956521+1400000], "chr5" => [11725524-500000, 11725524+500000]}

class Remove_Cen
	def self.tromere(cen, positions, chromosome)
		nocen = []
		positions.each do |pos|
			c_start, c_end = cen["chr#{chromosome}"]
			if pos < c_start || pos > c_end
				nocen << pos 
			end
		end
		return nocen
	end 
end 

hm_nocen = Remove_Cen.tromere(centromere, hm_list, chromosome)
ht_nocen = Remove_Cen.tromere(centromere, ht_list, chromosome)



File.open("Reads/#{location}/hm_nocen.txt", "w+") do |f|
  hm_nocen.each { |i| f.puts(i) } 
end


File.open("Reads/#{location}/ht_nocen.txt", "w+") do |f|
  ht_nocen.each { |i| f.puts(i) }
end

puts "creating the new files"
