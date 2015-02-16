#encoding: utf-8
class ModelGenome

	require 'rubygems'
	require 'bio-samtools'
	require 'bio'
	require 'rinruby'

	# Make a list that models SNP positions
	# Input 0: String of r code for the desired homozygous SNP distribution
	# Input 1: String of r code for the desired heterozygous SNP distribution
	# Output 0: List of homozygous SNPs
	# Output 1: List of heterozygous SNPs
	def self.get_snps(hm_code, ht_code)
		myr = RinRuby.new(echo = false)
		myr.eval hm_code
		myr.eval ht_code
		hm = myr.pull 'hm'
		ht = myr.pull 'ht'
		myr.quit
		return hm.map(&:abs).map(&:to_i).uniq, ht.map(&:abs).map(&:to_i).uniq # a few SNPs may be removed but doesn't affect distribution much
	end

	# Input: FASTA file
	# Output: Character array for all the bases in the first entry of the input FASTA
	def self.fasta_to_char_array (fasta)
		fasta_array = []
		x = 1
		Bio::FastaFormat.open(fasta).each do |i|
			fasta_array << i.seq
			# puts "f #{(x/fasta.length)*100}%"
			x+=1
		end
		return fasta_array[0].split(//)
	end

	# Input 0: Character array of genome sequence e.g. ['A','T', 'C'...]
	# Input 1: Size of fragments: Output fragments are of random length between the value of input 1 and double that value
	# Output: Fragmented genome: Character array split into chunks and put into super array e.g. [['A', 'T'], ['C', 'G']...]
	def self.get_frags (seq, size)
		frags = []
		rt, x = 0, 1
		while frags.flatten.length < seq.length
			frag_length = rand(size) + size # fragments of 10kb - 20kb
			frag = seq[rt..(rt + frag_length)]
			frags << frag
			rt = rt + frag_length
			x+=1
		end
		frags = frags[0..-2]
		if frags.flatten.length < seq.length
			frag = seq[frags.flatten.length..seq.length-1]
			frags << frag
		end
		return frags
	end

	# Input 0: Array of all the SNP positions in the genome
	# Input 1: Fragmented genome: Character array split into chunks and put into super array e.g. [['A', 'T'], ['C', 'G']...]
	# Output 0: Array of arrays where each sub array contains the SNP positions on one fragment, and each sub array is ordered in the same way as Input 1 in the super array.
	# Output 1: Array of arrays where each sub array contains the genomic positions of the SNPs for that fragment, and each sub array is ordered in the same way as Input 1 in the super array.
	def self.pos_each_frag (snp_pos, frags) # get the positions for snps on individual frags
		snp_pos.sort! # this is needed as the ht/hm SNPs need to be ordered together
		p_ranges = []
		frags.each do |i|
			p_ranges << i.length + p_ranges[-1].to_i # adding the fragment lengths to get the upper bounds of ranges of positions on the original seq.
		end
		first_pos = [] # then, to work out the first position of each fragment
		first_pos << 0
		first_pos << p_ranges[0..-2]
		first_pos = first_pos.flatten
		p = 0
		t = 0
		all_frags_pos = [] # the positions of snps on each fragment, array of arrays
		snp_pos_all = [] # the actual positions in the genome for each fragment
		p_ranges.each do |jj| # for each of the upper ranges (j) that the positions could be in
			each_frag_pos = []
			actual_pos = []
			while snp_pos[t].to_i < jj && snp_pos[t].nil? == false do # make the loop quit before trying to access index of snp_pos that doesn't exist
				each_frag_pos << snp_pos[t].to_i 	# add all of the positions < jj to a new array for each frag
				actual_pos << snp_pos[t].to_i
				t += 1
			end
			snp_pos_all << actual_pos 
			z = 0
			y = first_pos[p].to_i
			each_frag_pos.each do |e|
	      each_frag_pos[z] = e - y # taking away the value at the first position of each fragment to get the true/actual positions
	      z += 1
	    end
			p += 1
			all_frags_pos << each_frag_pos # add the positions for the frag to an array of the each_frag arrays
		end
		return all_frags_pos, snp_pos_all
	end

	###########
	#FASTA/VCF#
	###########

	# Input: Fragmented genome: Character array split into chunks and put into super array e.g. [['A', 'T'], ['C', 'G']...]
	# Output: Array of arrays, a sub array for each fragment that contains 2 strings: an identifier with the fragment's length, and the sequence of that fragment.
	def self.fasta_array (frags)
		frag_ids, id_and_length = [], []
		x = 1
		frag_strings = []
		frags.each do |i|
			id = "frag#{x}"
			id_and_length << ">#{id} Length = #{i.length}"
			frag_ids << id
			frag_strings << i.join.to_s # joining all the bases in one frag array to form a string for fasta
			x+=1
		end
		fastaformat_array = id_and_length.zip(frag_strings) #create the array, each element of which goes onto a new line in fasta
		return fastaformat_array
	end

	# Input 0: Fragmented genome: Character array split into chunks and put into super array e.g. [['A', 'T'], ['C', 'G']...]
	# Input 1: Array of arrays where each sub array contains the SNP positions on one fragment
	# Input 2: Array of arrays where each sub array contains the genomic positions of the SNPs for that fragment
	# Input 3: Array of the genomic SNP positions that are homozygous
	# Input 4: Array of the genomic SNP positions that are heterozygous
	# Output: Array of strings, each string is a line that will be written to the VCF file
	def self.vcf_array (frags, pos_on_frags, snp_pos_all, hm, ht)
		chrom, ref = [], []
		q = 1
		while q <= pos_on_frags.length
			if pos_on_frags[q-1].length != 0 #all of the fragments that contain at least one snp
				pos_on_frags[q-1].length.times do #*no. of snps
					chrom << "frag#{q}"
				end
			end
			pos_on_frags[q-1].each do |i|
				ref << frags[q-1][i] #what nucleotide is at these positions?
			end
			q+=1
		end
		alt = []
		it = 0
		ref.each do |base| 
			base == nil ? ref[it] ='T': case base
			when 'A' then alt << 'T'
			when 'T' then alt << 'A'
			when 'C' then alt << 'G'
			when 'G' then alt << 'C'
			when 'N' then alt << 'R'
			else alt << 'A'
			end
			it+=1
		end
		vcf_format = ['##fileformat=VCFv4.1', '##source=Fake', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'] 
		u = 0
		pos_on_frags.flatten.each do |i| #if we increment the pos_on_frags and the snp_pos_all together, we can tell whether each SNP from on_frags is het/homo
			if hm.include?(snp_pos_all.flatten[u])
				x = '1.0'
			elsif ht.include?(snp_pos_all.flatten[u])
				x = '0.5'
			else
				x = 'WRONG'
			end
			line = "#{chrom[u]}	#{i}	.	#{ref[u]}	#{alt[u]}	100	PASS	AF=#{x}"
			vcf_format << line
			u+=1
		end
		return vcf_format
	end
end