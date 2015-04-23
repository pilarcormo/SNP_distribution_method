
#encoding: utf-8
require_relative 'lib/write_it'

###get rid of centromeres: 500 kb at each side of the centromere
chromosome = ARGV[0] 
location = ARGV[1]

hm_list = WriteIt.file_to_ints_array("Reads/#{location}/hm.txt") # Get SNP distributions
ht_list = WriteIt.file_to_ints_array("Reads/#{location}/ht.txt")

cen_1 = Array(15086545 - 500000..15086545+500000)
cen_2 = Array(3608429- 500000..3608429+ 500000)
cen_3 = Array(14209452- 500000..14209452+500000)
cen_4 = Array(3956521- 500000..3956521+500000)
cen_5 = Array(11725524-500000..11725524+500000)

class Remove_Cen
	def self.tromere(cen, hm_list)
		hm_nocen = []
		hm_list.each do |pos|
			if cen.include?(pos)
			else 
				hm_nocen << pos 
			end
		end
		return hm_nocen
	end 
end 

if chromosome.to_i == 1
	hm_nocen = Remove_Cen.tromere(cen_1, hm_list)
	ht_nocen = Remove_Cen.tromere(cen_1, ht_list)
end 
if chromosome.to_i == 2
	hm_nocen = Remove_Cen.tromere(cen_2, hm_list)
	ht_nocen = Remove_Cen.tromere(cen_2, ht_list)
end 
if chromosome.to_i == 3
	hm_nocen = Remove_Cen.tromere(cen_3, hm_list)
	ht_nocen = Remove_Cen.tromere(cen_3, ht_list)
end 
if chromosome.to_i == 4
	hm_nocen = Remove_Cen.tromere(cen_4, hm_list)
	ht_nocen = Remove_Cen.tromere(cen_4, ht_list)
end 
if chromosome.to_i == 5
	hm_nocen = Remove_Cen.tromere(cen_5, hm_list)
	ht_nocen = Remove_Cen.tromere(cen_5, ht_list)
end 


File.open("Reads/#{location}/hm_nocen.txt", "w+") do |f|
  hm_nocen.each { |element| f.puts(element) }
end
File.open("Reads/#{location}/ht_nocen.txt", "w+") do |f|
  ht_nocen.each { |element| f.puts(element) }
end


