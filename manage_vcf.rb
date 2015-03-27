#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require 'pp'

file = ARGV[0]
chromosome = ARGV[1].to_i
vcf = ARGV[2]

vcf_file = "#{file}/#{vcf}.vcf"
location = "#{file}/chromosome#{chromosome}"

new_vcf, vaf = [], []
Dir.mkdir(location)

File.open(vcf_file, 'r').each do |line|
	next if line =~ /^#/
    v = Bio::DB::Vcf.new(line)
    a = line.split("\t")
    first = a.first
    if first == "#{chromosome}"
    	new_vcf << line 
    end
end 

File.open("#{location}/vcf_file_#{chromosome}.vcf", "w+") do |f|
  new_vcf.each { |element| f.puts(element) }
end

vcf_file_chr = "#{location}/vcf_file_#{chromosome}.vcf"

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

WriteIt::write_txt("#{location}/hm", hm) 
WriteIt::write_txt("#{location}/ht", ht)