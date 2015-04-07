
#encoding: utf-8
class Vcf 
	def self.type_per_pos (vcfs_info, vcfs_pos)
	    snps = {}
	    x = 0
	    vcfs_info.each do |hash|
	    	hash.each do |type, number|
	    		if number == "1" 
	    			snps.store(vcfs_pos[x], type) 		
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
	    return snps, hm, ht
	end  
    def self.open_vcf(vcf_file, chromosome)
	    vcfs_chrom, vcfs_pos, vcfs_info = [], [], [] 
	    new_vcf = []
	    File.open(vcf_file, 'r').each do |line|
	    	next if line =~ /^#/
	        v = Bio::DB::Vcf.new(line)
	        vcfs_chrom << v.chrom
	        vcfs_pos << v.pos
	        vcfs_info << v.info 
	        # snp_type = v.info["HOM"].to_f
	        a = line.split("\t")
	        if v.chrom == "#{chromosome}"
	        	new_vcf << line 
	        end
	    end 
		return new_vcf, vcfs_chrom, vcfs_pos, vcfs_info
	end 
end 
