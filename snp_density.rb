require_relative "lib/write_it"
require 'pp'
require 'csv'
number = ARGV[0]
folder = ARGV[1]
pre_filter = ARGV[2]
post_filter = ARGV[3]
name = ARGV[4]

hm = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/hm.txt")
ht = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/ht.txt")

pp hm.length 
pp ht.length

hm_fil = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/#{post_filter}/hm.txt")
ht_fil = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/#{post_filter}/ht.txt")

pp hm_fil.length 
pp ht_fil.length


hm_noc = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/#{post_filter}/hm_nocen.txt")
ht_noc = WriteIt.file_to_array("Reads/#{folder}/#{pre_filter}/#{post_filter}/ht_nocen.txt")

pp hm_noc.length 
pp ht_noc.length


if File.exist?('density.csv') == true 
  CSV.open("density.csv", "a+") do |csv|
    csv << ["#{name}", number, hm.length, hm_fil.length, hm_noc.length]
  end 
else 
  CSV.open("density.csv", "wb") do |csv|
    csv << ["study", "chromosome", "pre-filtering", "parental", "centromere"]
  end 
  CSV.open("density.csv", "a+") do |csv|
    csv << ["#{name}", number, hm.length, hm_fil.length, hm_noc.length]
  end 
end 
