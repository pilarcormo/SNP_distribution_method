#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/vcf'
require 'pp'

file = ARGV[0]
vcf = ARGV[1]
chromosome = ARGV[2].to_i
filtering = ARGV[3]
parent = ARGV[4]

location = "#{file}/#{vcf}_chromosome#{chromosome}"


#######

case filtering
when 'background_vcf' 
    puts "Opening the vcf file"

    Dir.mkdir("#{location}")
    vcf, vcfs_chrom, vcfs_pos, vcfs_info = Vcf.open_vcf("#{file}/#{vcf}.vcf", chromosome)
    
    vcf_file_chr = "#{location}/chromosome#{chromosome}.vcf"
    puts "Creating VCF file for chromosome #{chromosome}"
    
    File.open(vcf_file_chr, "w+") do |f|
      vcf.each { |element| f.puts(element) }
    end

    chr_vcf, vcfs_chrom_p, vcfs_pos_p, vcfs_info_p = Vcf.open_vcf(vcf_file_chr, chromosome)
    snps_p, hm_p, ht_p = Vcf.type_per_pos(vcfs_info_p, vcfs_pos_p)

    WriteIt::write_txt("#{location}/hm", hm_p) 
    WriteIt::write_txt("#{location}/ht", ht_p)

when 'filter_vcf' 
    puts "Opening the child vcf file"
    vcf_file_chr = "#{location}/chromosome#{chromosome}.vcf"
    chr_vcf, vcfs_chrom_c, vcfs_pos_c, vcfs_info_c = Vcf.open_vcf(vcf_file_chr, chromosome)

    snps_c, hm_c, ht_c = Vcf.type_per_pos(vcfs_info_c, vcfs_pos_c)
    puts "Opening the parental vcf file"
    parental_chr_vcf, vcfs_chrom_p, vcfs_pos_p, vcfs_info_p = Vcf.open_vcf("#{parent}", chromosome)

    snps_p, hm_p, ht_p = Vcf.type_per_pos(vcfs_info_p, vcfs_pos_p)

    short_child_chr_vcf = Vcf.filtering(vcfs_pos_c, snps_p, snps_c, chr_vcf)
    puts "Writing the interesting SNPs"

    int = "interesting_#{chromosome}"
    Dir.mkdir(int)

    File.open("#{int}/chromosome#{chromosome}_interesting.vcf", "w+") do |f|
      short_child_chr_vcf.each { |element| f.puts(element) }
    end

    interesting = "#{int}/chromosome#{chromosome}_interesting.vcf"

    int_chr_vcf, vcfs_chrom_i, vcfs_pos_i, vcfs_info_i = Vcf.open_vcf(interesting, chromosome)
    snps_i, hm_i, ht_i = Vcf.type_per_pos(vcfs_info_i, vcfs_pos_i)

    WriteIt::write_txt("#{int}/hm", hm_i) 
    WriteIt::write_txt("#{int}/ht", ht_i)
else 
    puts "lalalala"
end 



######





