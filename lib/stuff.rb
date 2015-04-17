
#encoding: utf-8
require_relative 'reform_ratio'
require 'pp'

##Open the vcf file and create lists of heterozygous and homozygous SNPs

class Stuff
	##Input: vcf file
	##Ouput: lists of hm and ht SNPS
	def self.snps_in_vcf(vcf_file, ht_cutoff=0.5, hm_cutoff=1.0)
		vcfs_chrom, vcfs_pos, vcfs_info, hm, ht = [], [], [], [], []
		frag_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
		File.open(vcf_file, "r").each do |line| # get array of vcf lines, you can call a method on one line
			next if line =~ /^#/
			v = Bio::DB::Vcf.new(line)
			vcfs_chrom << v.chrom
			vcfs_pos << v.pos
			vcfs_info << v.info # so this will be an array of hashes of strings
			allele_freq = v.info["AF"].to_f
			if allele_freq == ht_cutoff
		        if frag_pos[:het].has_key?(v.chrom)
		        	frag_pos[:het][v.chrom] << v.pos
		        else
		        	frag_pos[:het][v.chrom] = []
		        	frag_pos[:het][v.chrom] << v.pos
		        end
		        ht << v.chrom
		    elsif allele_freq == hm_cutoff
		       	if  frag_pos[:hom].has_key?(v.chrom)
		        	frag_pos[:hom][v.chrom] << v.pos
		        else
		        	frag_pos[:hom][v.chrom] = []
		        	frag_pos[:hom][v.chrom] << v.pos
		        end
		        hm << v.chrom 
		    end
		end 
		num_snps_frag_hash = Hash.new(0)
		vcfs_chrom.each {|v| num_snps_frag_hash[v] +=1 } # we have the number of snps on each frag, by counting the repeats of each frag in the vcf
		# # the frag_id(.chrom) is the key, the number of snps for that frag is the value. putting the number of snps for each frag into hash
		snp_data = vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info
		return snp_data, hm, ht
	end

	def self.safe_invert(hash)
    	hash.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
 	end

	def self.dic_id_pos(h, snp_list)
		dic_pos = {}
	  	x = 0 
	  	Array(0..snp_list.length - 1).each do |o|
	    	dic_pos.store(snp_list[x], h[x])
	   	 	x += 1 
	 	end
	 	dic_pos = Stuff::safe_invert(dic_pos)
	  	return dic_pos
	end 



	##Input: Lists of hm and ht SNPs
	##Output: dictionaries with the id of the fragment as key and the absolute number of SNPs as value
	def self.create_hash_number(array)
		hash1 = {}
		array.uniq.each { |elem| hash1.store("#{elem}", "#{array.count(elem).to_i}") }
		return hash1
	end
	##Input 1: Array of fragment ids.
	##Input 2: Hash of hm SNPs
	##Input 3: Hash of ht SNPs
	##Assign the number of SNPs to each fragment.
	##If a fragment does not have SNPs, the value assigned will be 0.
	##Output1: New hash with the number of SNPs assigned to the unordered fragments
	##Output2: Array with the number of SNPs 
  def self.define_snps(ids, dic_snps)
		shuf = {}
		snps = []
		ids.each { |frag|
	      if dic_snps.has_key?(frag)
	        shuf.store(frag, dic_snps[frag].to_f)
	      else
	        shuf.store(frag, 0)
	      end
	  	}
		shuf.each { |id, snp| snps << snp }
		return shuf, snps
	end

	##Inputs: hashes with IDs as keys and the SNP density as value
	##Divide absolute number of SNPs by the length of the given fragment. 
	##Output: hashes with IDs as keys and the normalised SNP density as value
	def self.normalise_by_length(lengths, shuf)
		shuf_norm = {}
		l = shuf.length
		Array(0..l-1).each { |x|
	      snp_norm = shuf[shuf.keys[x]].to_f/lengths[x].to_f
	      shuf_norm.store(shuf.keys[x], snp_norm)
    	}
		shuf_norm = Stuff::safe_invert(shuf_norm)
		return shuf_norm
	end

	##Input1: permutation array after SDM
	##Input2: fasta file converted to array
	##Input3: list of fragment ids from the shuffled fasta file
	##Output: permutation of fragments after SDM with the data (lengths, etc) obtained from the original fasta file

	def self.create_perm_fasta(perm, frags, ids)
		defs, data, defs_p, data_p, fasta_perm = [], [], [], [], []
		frags.each do |i|
			defs << i.definition
			data << i.data
		end
		perm.each { |frag|
	      index_frag = ids.index(frag).to_i
	      defs_p << defs[index_frag]
	      data_p << data[index_frag]
    	}
    ###Create fasta array with the information above
		x = 0
		Array(0..perm.length-1).each do |i|
			fasta_perm << defs_p[x]
			fasta_perm << data_p[x]
			x += 1
		end
		fasta_perm.each do |line|
			if line.start_with?("\n")
				line[0] = ''
			else
				line.insert(0, ">")
			end
		end
		return fasta_perm
	end

	def self.positions_by_fragment(dic, snp_list)
		dic.each do |id, number|
			dic.store(id, snp_list[0..(number.to_i-1)])
			(number.to_i).times do
		    	snp_list.delete_at(0)
			end
		end
		dic.delete_if { |id, number|  number.empty?}
		return dic
	end

	def self.important_ratios(snps_hm, snps_ht, ids, id_len, threshold, adjust) 
		x = 0
		dic_ratios, ratios, ids_s, ratios_by_len = {}, [], [], []
		dic_ratios_by_length = {}
		snps_hm.length.times do
			ratio = (snps_hm[x]+adjust.to_f)/(snps_ht[x]+adjust.to_f)
			dic_ratios.store(ids[x], ratio.to_f) 
			x += 1
		end
		if threshold.to_i > 0 
			dic_ratios.delete_if { |id, ratio|  ratio <= threshold.to_f}
		end 
		ratios << dic_ratios.values
		ids_s << dic_ratios.keys
		ratios.flatten!
		x = 0 
		id_len.each do |id, len|
			if dic_ratios.has_key?(id)
				ratio_by_len = ratios[x].to_f/len.to_i
				ratios_by_len << ratio_by_len
				dic_ratios_by_length.store(id, ratio_by_len)
			end 
			x += 1
		end 
		return dic_ratios, ratios, ids_s, ratios_by_len, dic_ratios_by_length
	end

	def self.important_ids(ids_short, ids)
		shuf_short_ids = []
		ids_short.flatten!
		ids.each do |frag|
			if ids_short.include?(frag)
		    	shuf_short_ids << frag
		  	end 
		end 
		shuf_short_ids.flatten!
		return shuf_short_ids
	end 
	def self.important_pos(ids_short, pos)
		sh = []
		pos.each do |frag, positions|
			if ids_short.include?(frag)
		  	else 
		    	pos.delete(frag)
		  	end 
		end 
		sh = pos.values
		sh.flatten!
		return sh
	end 

	#Input1 location for the csv file
	#Input2 hash with the id and positions for the hm SNPs
	#Input3 hash with the id and ratio for each fragment
	def self.csv_pos_ratio(csv, pos, ratios)
		CSV.open(csv, "wb") do |csv|
		  csv << ["Position", "Ratio"]
		end
		short = {}
		short = pos
		short.each do |id, array|
		  if ratios.has_key?(id)
		  else 
		    short.delete(id)
		  end
		end
		short.each do |id, array|
			array.each do |elem| 
		    	CSV.open(csv, "ab") do |csv|
		      		csv << [elem, ratios[id]] 
		      	end 
		    end 
		end 
	end 
end