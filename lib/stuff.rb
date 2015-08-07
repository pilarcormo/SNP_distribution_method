

#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pp'

##Open the vcf file and create lists of heterozygous and homozygous SNPs

class Stuff
	# Input: Array of Bio::FastaFormat entries
	# Output 0: Array of identifiers
	# Output 1: Array of lengths (integers)
	def self.fasta_id_n_lengths(fasta)
		ids, lengths = [], []
		id_len = {}
		fasta.each do |i|
			ids << i.entry_id
			lengths << i.length
			id_len.store(i.entry_id, i.length)
		end
		return ids, lengths, id_len
	end

	# Input: FASTA file
	# Output: Array of Bio::FastaFormat entries
	def self.fasta_array(fasta_file)
		fasta = [] # we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
		Bio::FastaFormat.open(fasta_file).each do |i| # get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
			fasta << i
		end
		return fasta
	end

	# Input: FASTA file
	# Output: Integer of the length of the genome
	def self.genome_length(fasta_file)
		lengths = []
		fasta_array(fasta_file).each do |frag|
			lengths << frag.length
		end
		return lengths.inject(:+)
	end
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
		return snp_data, hm, ht, frag_pos
	end

	def self.safe_invert(hash)
    	hash.each_with_object( {} ) { |(key, value), out| ( out[value] ||= [] ) << key }
 	end

	def self.dic_id_pos(h_ids, snp_list)
		dic_pos = {}
	  	x = 0 
	  	Array(0..snp_list.length - 1).each do |o|
	  		if dic_pos.has_key?(h_ids[x])
	  			dic_pos[h_ids[x]] << snp_list[x]
	  		else 
	  			dic_pos[h_ids[x]] = []
	  			dic_pos[h_ids[x]] << snp_list[x]
	  		end 
	    	x += 1 
	 	end
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
		dic_snps_num = {}
		snps = []
		ids.each { |frag|
	      if dic_snps.has_key?(frag)
	        dic_snps_num.store(frag, dic_snps[frag].to_f)
	      else
	        dic_snps_num.store(frag, 0)
	      end
	  	}
		dic_snps_num.each { |id, snp| snps << snp }
		return dic_snps_num, snps
	end

	def self.define_global_pos(ids, dic_pos, lengths)
		global_array, lens = [], []
	  	or_pos, temporal_pos, dic_global_pos, idlen_s = {}, {}, {}, {}

	  	pp ids, dic_pos, lengths
		# ids.each { |frag|
		#     if dic_pos.has_key?(frag)
		#       or_pos.store(frag, dic_pos[frag])
		#     end
	 #  	}
	 #  	keys = or_pos.keys
	  	# temporal_pos.store(keys[0], or_pos[keys[0]])
	  	# or_pos.delete_if { |frag, pos|  temporal_pos.has_key?(frag)}
	  	# id_length.each do |frag, length|
	  	# 	idlen_s.store(frag, length)
	  	# 	lens << length
	  	# end 
	  	x = 1 
	  	lengths.length.times do |i|
	  		lengths[x] = lengths[x-1].to_i + lengths[x].to_i 
	  		x += 1
	  	end 

	  	x = 0 
    	dic_pos.each do |frag, array|
	  		array.each do |pos|
	  			pos2 = (pos.to_i+lengths[x].to_i)
	  			global_array << pos2
	  		end  
	  		x += 1
	  		dic_global_pos.store(frag, global_array)
  			global_array = []
  		end 
	  	#dic_global_pos = temporal_pos.merge(dic_global_pos)
	  	all_global_positions = dic_global_pos.values
	  	all_global_positions.flatten!
		return dic_global_pos, all_global_positions
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

end