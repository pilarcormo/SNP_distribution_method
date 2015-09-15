#encoding: utf-8
class ReformRatio
	require 'bio'
	require 'bio-samtools'
	# Input: Array of Bio::FastaFormat entries
	# Output 0: Array of identifiers
	# Output 1: Array of lengths (integers)
	def self.fasta_id_n_lengths(fasta)
		ids, lengths = [], []
		fasta.each do |i|
			ids << i.entry_id
			lengths << i.length
		end
		return ids, lengths
	end

	# Input: VCF file
	# Output 0: Array of VCF chrom field (fragment identifiers)
	# Output 1: Array of VCF pos field (positions of snps on fragments)
	# Output 2: Hash with fragment id keys and the corresponding number of snps as an integer values
	# Output 3: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency)
	def self.get_snp_data(vcf_file)
		vcfs_chrom, vcfs_pos, vcfs_info = [], [], []
		File.open(vcf_file, "r").each do |line| # get array of vcf lines, you can call a method on one line
			next if line =~ /^#/
			v = Bio::DB::Vcf.new(line)
			vcfs_chrom << v.chrom
			vcfs_pos << v.pos
			vcfs_info << v.info # so this will be an array of hashes of strings
		end
		num_snps_frag_hash = Hash.new(0)
		vcfs_chrom.each {|v| num_snps_frag_hash[v] +=1 } # we have the number of snps on each frag, by counting the repeats of each frag in the vcf
		# the frag_id(.chrom) is the key, the number of snps for that frag is the value. putting the number of snps for each frag into hash
		return vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info
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

	# Input 0: Hash with fragment id keys and the corresponding number of snps as an integer values
	# Input 1: Array of Bio::FastaFormat entries
	# Output: Array of the number of snps per fragment, in the same order as the input 1 fasta array
	def self.snps_per_fasta_frag(snps_per_vcf_frag_hash, fasta_array)
		snps_per_frag_fasta_order = [] #use the id to identify the number of snps for that frag using the keys of snps_hash
		fasta_array.each do |frag|
			snps_per_frag_fasta_order << snps_per_vcf_frag_hash[frag.entry_id] #gives 0 for keys that don't exist = good, because the frags with 0 density would otherwise be missing
		end #now we have an array with the number of snps per frag in the same order as the fasta array
		return snps_per_frag_fasta_order
	end

	# Input 0: Array of Bio::FastaFormat entries
	# Input 1: Array of VCF chrom field (fragment identifiers)
	# Input 2: Array of VCF pos field (positions of snps on fragments)
	# Input 3: Array of the number of snps per fragment, in the same order as the input 0 fasta array
	# Input 4: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency, single key/value)
	# Output 0: The snp positions for each frag, in an array of arrays (each sub array contains the snp positions for one frag, and the sub arrays are ordered according to the order of the input 1 fasta)
	# Output 1: The info hashes (of each snp, single key/value) for each frag, in an array of arrays (each sub array contains the info hashes for one frag, and the sub arrays are ordered according to the order of the input 1 fasta)
	def self.get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag, vcfs_info)
		pos, info = [], []
		n = 0
		fasta.each do |frag|
			x = 0
			each_fr_pos, each_fr_info = [], []
			snps_per_frag[n].times do |i|
				while frag.entry_id != vcfs_chrom[x]
					x+=1
				end
				each_fr_pos << vcfs_pos[x]
				each_fr_info << vcfs_info[x]
				x+=1
			end
			pos << each_fr_pos #this gives empty arrays for frags with out snps, and a list of the positions of those with
			info << each_fr_info
			n+=1
		end
		return pos, info
	end

	# Input 0: The snp positions for each frag, in an array of arrays (for a given fragment permutation)
	# Input 1: Array of fragment lengths (integers) (for the same fragment permutation as Input 0)
	# Output: Array of the snp positions in the genome (assuming the genome is ordered according to the fragment permutation in Input 0)
	def self.total_pos(pos, fasta_lengths)  # both args in same order as fasta = good. this all works!
		totals = []                    
		x = 0						   
		pos.each do |frag|
			if x == 0
				totals << frag
				x+=1
			else
				tot_frag, lengths = [], []
				fasta_lengths[0..x-1].each do |l|
					lengths << l
				end
				so_far = lengths.inject(:+) # this is the sum of the lengths for all the frags so far
				frag.each do |i|
					tot_frag << so_far + i
				end
				totals << tot_frag
				x+=1
			end
		end
		return totals.flatten
	end

	# Input 0: Array of the snp positions in the genome
	# Input 1: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency)
	# Output 0: Array of all the heterozygous snp positions
	# Output 1: Array of all the homozygous snp positions
	def self.het_hom(actual_pos, vcfs_info) #actual_pos in same order as fasta perm. now so is info
		het, hom = [], []
		x = 0
		actual_pos.each do |snp|
			if vcfs_info.flatten[x]['AF'] == '1.0' # homozygous SNPs have AF= 1.0, we can change this to a range for real data
				hom << snp
			elsif vcfs_info.flatten[x]['AF'] == '0.5'
				het << snp
			end
			x+=1
		end
		return het, hom
	end

	# Input 0: Array of FASTA format objects (contig permutation)
	# Input 1: SNP data output from get_snps method
	# Output 0: Array of the heterozyogus SNP positions for the genome according to the Input 0 permutation
	# Output 1: Array of the homozygous SNP positions for the genome according to the Input 0 permutation
	def self.perm_pos(fasta, snp_data)
		snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) # array of no. of snps per frag in same order as fasta
		pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) # get snp positions for each frag in array of arrays
		actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
		het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])
		return het_snps, hom_snps
	end
end