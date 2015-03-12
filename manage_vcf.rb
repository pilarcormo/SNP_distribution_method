#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require 'pp'

vcf_file = "BCF2.vcf"
# vcf_file = "Varscan.vcf"
chromosome = ARGV[0].to_i

m = [] 
Dir.mkdir("BCF2_chromosome#{chromosome}_nocentromere")

File.open(vcf_file, 'r').each do |line|
	next if line =~ /^#/
    v = Bio::DB::Vcf.new(line)
    a = line.split("\t")
    first = a.first
    if first == "#{chromosome}"
    	m << line 
    end
end 

File.open("BCF2_chromosome#{chromosome}_nocentromere/vcf_file_#{chromosome}.vcf", "w+") do |f|
  m.each { |element| f.puts(element) }
end

vcf_file_chr = "BCF2_chromosome#{chromosome}_nocentromere/vcf_file_#{chromosome}.vcf"

vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info = ReformRatio.get_snp_data(vcf_file_chr)

snps = {}
x = 0

vcfs_info.each do |hash|
	hash.each do |thing, number|
		if number == "1" 
			snps.store(vcfs_pos[x], thing) 		
			x += 1		
		end
	end
end 
 	
hm, ht = [], [] 
snps.each do |pos, type|
	if type == "HET"
		ht << pos
	elsif type == "HOM"
		hm << pos
	end 
end 	

hm2 = hm 

hm_rm, ht_rm = [], []


hm.each do |item|
	item = item.to_f/100000
	if item <= 151 && item >= 149 
		hm.delete(item)
		hm_rm << item 
	end
end

pp hm_rm.length
pp hm.length

# hm.each do |item|
# 	if item.to_i <= 15100000 && item.to_i >= 14900000
# 		hm.delete(item)
# 	end
# ht.each do |item|
# 	if item.to_i <= 15100000 && item.to_i >= 14900000
# 		ht.delete(item)
# 	end

hm_rm.each do |y|
	pp y 
	y = y*100000
	hm2.delete(y)
end

pp hm2.length




# WriteIt::write_txt("BCF2_chromosome#{chromosome}_nocentromere/hm", hm) 
# WriteIt::write_txt("BCF2_chromosome#{chromosome}_nocentromere/ht", ht)


