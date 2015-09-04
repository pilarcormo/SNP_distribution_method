#encoding: utf-8
require_relative 'lib/write_it'
require_relative 'lib/vcf'
require 'pp'

if ARGV.empty?
    puts "Please specify (1) a location folder,  (2) a vcf file, (3) the chromosome we want to 
    analyse (1, 2, 3, 4, 5..) and (4) a method for managing the vcf file: 
    use either 'cutting_vcf' or 'filter_vcf'. If the second option is used we also need (5) a parental vcf file"
else 
    file = ARGV[0] #location folder where the vcf file is and where the output vcf files will be created.
    vcf = ARGV[1] #starting vcf file 
    folder = ARGV[2] #output folder
    chromosome = ARGV[3] #chromosome we want to analyse
    filtering = ARGV[4] #either cutting_vcf (to divide the original vcf file by chromosome in individual vcf files and create hm.txt and ht.txt) 
    #or filter_vcf (use the first file as child vcf and the parent provided after as filter to eliminate brackground SNPs)       
    parent = ARGV[5]
end 
location = "#{file}"
output = "#{file}/#{folder}"

case filtering

when 'cutting_vcf'
    Dir.mkdir(output)

    puts "Opening the vcf file"
    vcf, vcfs_chrom, vcfs_pos, vcfs_info = Vcf.open_vcf("#{file}/#{vcf}", chromosome)

    vcf_file_chr = "#{output}/chromosome#{chromosome}.vcf"
    puts "Creating VCF file for chromosome #{chromosome}"
    
    File.open(vcf_file_chr, "w+") do |f|
      vcf.each { |element| f.puts(element) }
    end
    puts "Obtaining the lists of hm and ht SNPs and saving them in text files"
    chr_vcf, vcfs_chrom_p, vcfs_pos_p, vcfs_info_p = Vcf.open_vcf(vcf_file_chr, chromosome)

    snps_p, hm_p, ht_p = Vcf.type_per_pos(vcfs_info_p, vcfs_pos_p)

    WriteIt::write_txt("#{output}/hm", hm_p) 
    WriteIt::write_txt("#{output}/ht", ht_p)

when 'filter_vcf' 
    puts "Opening the child vcf file"
    vcf_file_chr = "#{location}/#{vcf}"
    chr_vcf, vcfs_chrom_c, vcfs_pos_c, vcfs_info_c = Vcf.open_vcf(vcf_file_chr, chromosome)

    snps_c, hm_c, ht_c = Vcf.type_per_pos(vcfs_info_c, vcfs_pos_c)
    puts "Opening the parental vcf file"
    parental_chr_vcf, vcfs_chrom_p, vcfs_pos_p, vcfs_info_p = Vcf.open_vcf("#{file}/#{parent}", chromosome)

    snps_p, hm_p, ht_p = Vcf.type_per_pos(vcfs_info_p, vcfs_pos_p)

    short_child_chr_vcf = Vcf.filtering(vcfs_pos_c, snps_p, snps_c, chr_vcf)

    puts "Writing the interesting SNPs"

    int = "#{location}/#{folder}"
    Dir.mkdir(int)

    File.open("#{int}/chromosome#{chromosome}.vcf", "w+") do |f|
      short_child_chr_vcf.each { |element| f.puts(element) }
    end

    interesting = "#{int}/chromosome#{chromosome}.vcf"

    int_chr_vcf, vcfs_chrom_i, vcfs_pos_i, vcfs_info_i = Vcf.open_vcf(interesting, chromosome)
    snps_i, hm_i, ht_i = Vcf.type_per_pos(vcfs_info_i, vcfs_pos_i)

    WriteIt::write_txt("#{int}/hm", hm_i) 
    WriteIt::write_txt("#{int}/ht", ht_i)
else 
    puts "A method for managing the vcf file needs to be provided. Use either 'cutting_vcf' or 'filter_vcf'"
end 


