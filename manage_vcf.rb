#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require 'pp'

vcf_file = "BCF2.vcf"
# vcf_file = "Varscan.vcf"
chromosome = ARGV[0]

m = [] 
Dir.mkdir("BCF2_chromosome#{chromosome}")

File.open(vcf_file, 'r').each do |line|
	next if line =~ /^#/
    v = Bio::DB::Vcf.new(line)
    a = line.split("\t")
    first = a.first
    if first == "#{chromosome}"
    	m << line 
    end
end 

File.open("BCF2_chromosome#{chromosome}/vcf_file_#{chromosome}.vcf", "w+") do |f|
  m.each { |element| f.puts(element) }
end

vcf_file_chr = "BCF2_chromosome#{chromosome}/vcf_file_#{chromosome}.vcf"

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

WriteIt::write_txt("BCF2_chromosome#{chromosome}/hm", hm) 
WriteIt::write_txt("BCF2_chromosome#{chromosome}/ht", ht)


