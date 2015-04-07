#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/vcf'
require 'pp'

file = ARGV[0]
chromosome = ARGV[1].to_i
vcf = ARGV[2]

parent_vcf = "#{file}/#{vcf}.vcf"
location = "#{file}/chromosome#{chromosome}"
Dir.mkdir(location)

child_vcf = "BCF2/BCF2_chromosome#{chromosome}/vcf_file_#{chromosome}.vcf"


parent_chr_vcf, vcfs_chrom_p, vcfs_pos_p, vcfs_info_p = Vcf.open_vcf(parent_vcf, chromosome)

snps_p, hm_p, ht_p = Vcf.type_per_pos(vcfs_info_p, vcfs_pos_p)

vcf_file_chr = "#{location}/vcf_file_#{chromosome}.vcf"

#write parental vcf file
File.open("#{location}/vcf_file_#{chromosome}.vcf", "w+") do |f|
  parent_chr_vcf.each { |element| f.puts(element) }
end

WriteIt::write_txt("#{location}/hm", hm_p) 
WriteIt::write_txt("#{location}/ht", ht_p)


child_chr_vcf, vcfs_chrom_c, vcfs_pos_c, vcfs_info_c = Vcf.open_vcf(child_vcf, chromosome)

snps_c, hm_c, ht_c = Vcf.type_per_pos(vcfs_info_c, vcfs_pos_c)


short_vcfs_pos_c = vcfs_pos_c
short_vcfs_pos_c.flatten!
snps_p.each do |pos, type|
    if snps_c.has_key?(pos)
        snps_c.delete(pos) 
        short_vcfs_pos_c.delete(pos)
    end 
end 

short_child_chr_vcf = []
child_chr_vcf.each do |line|
    puts  line 
    position = line.split("\t")[1].to_i
    puts position 
    if short_vcfs_pos_c.include?(position) 
        short_child_chr_vcf << line
    end 
end 

pp short_child_chr_vcf

int = "BCF2/Interesting_#{chromosome}"
Dir.mkdir(int)

File.open("#{int}/chromosome#{chromosome}_interesting.vcf", "w+") do |f|
  short_child_chr_vcf.each { |element| f.puts(element) }
end

interesting = "#{int}/chromosome#{chromosome}_interesting.vcf"

int_chr_vcf, vcfs_chrom_i, vcfs_pos_i, vcfs_info_i = Vcf.open_vcf(interesting, chromosome)
snps_i, hm_i, ht_i = Vcf.type_per_pos(vcfs_info_i, vcfs_pos_i)

WriteIt::write_txt("#{int}/hm", hm_i) 
WriteIt::write_txt("#{int}/ht", ht_i)

