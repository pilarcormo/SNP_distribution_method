class FitnessScore
	require 'rinruby'
	require_relative 'snp_dist'

	### Count/ratio method ################################################################################
		# Input 0: Array of SNP positions
		# Input 1: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 2: The length of the genome
		# Output: Array of number of SNPs in each genome division
		def self.count(snp_pos, div, genome_length)
			myr = RinRuby.new(echo = false)
			myr.assign 'snp_pos', snp_pos
			myr.assign 'div', div.to_i
			myr.assign 'l', genome_length
			myr.eval 'breaks <- c(0)
			for(i in 1:div){
			  breaks <- c(breaks,(l/div)*i)
			}'
			myr.eval 'counts <- hist(snp_pos, breaks=breaks, plot=FALSE)$counts'
			counts = myr.pull 'counts'
			myr.quit
			return counts
		end

		# Input 0: Array of homozygous SNP positions
		# Input 1: Array of heterozygous SNP positions
		# Input 2: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 3: The length of the genome
		# Output: Array of ratios (floats) of homozygous to heterozygous SNPs for each division of the genome permutation
		def self.ratio(hm, ht, div, genome_length)
			hm_count = FitnessScore::count(hm, div, genome_length)
			ht_count = FitnessScore::count(ht, div, genome_length)
			# pp hm_count
			# pp ht_count
			x = 0
			ratios = []
			div.times do
				count_ratio = ((hm_count[x]+1).to_f / (ht_count[x]+1).to_f) # a measure of ratio
				ratios << count_ratio
				x+=1
			end
			return ratios
		end

		# Input 0: Array of homozygous SNP positions
		# Input 1: Array of heterozygous SNP positions
		# Input 2: Number of breaks (divisions) in the genome to count the number of SNPs in
		# Input 3: The length of the genome
		# Input 4: Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division of the genome
		# Output: Float between 0.0 and 1.0 where closely matching inputs are closer to 1.0 (pearson correlation)
		def self.count_ratio(hm, ht, div, genome_length, expected_ratios)
			ratios = ratio(hm, ht, div, genome_length)
			myr = RinRuby.new(echo = false)
			myr.assign 'x', expected_ratios
			myr.assign 'y', ratios
			myr.eval 'score <- abs(cor(x,y))'
			fitness_score = myr.pull 'score'
			myr.quit
			return fitness_score
		end
	####################################################################################################

	# Input: Array of homozygous snp positions
	# Output: Integer of the total distance in bases, between adjacent SNPs
	def self.snp_distance(hm)
		score = 0
		hm.each_cons(2).map { |a,b| score+=(b-a) }
		return score
	end

	# Input: Array of homozygous snp positions
	# Output: Float of the maximum kernel density value of the homozygous SNP distribution
	def self.max_density(hm)
		myr = RinRuby.new(echo=false)
		myr.hm = hm
		myr.eval 'score <- max(density(hm)$y)'
		score = myr.pull 'score'
		myr.quit
		return score
	end

	# Input 0: Array of homozygous snp positions
	# Input 1: Array of heterozygous snp positions
	# Output: Float of the maximum ratio of homozygous to heterozygous kernel density
	def self.max_ratio(hm, ht)
		myr = RinRuby.new(echo=false)
		myr.hm, myr.ht = hm, ht
		myr.eval 'score <- max(density(hm)$y/density(ht)$y)'
		score = myr.pull 'score'
		myr.quit
		return score
	end

	# Input 0: Array of homozygous snp positions
	# Input 1: Array of heterozygous snp positions
	# Input 2: Number of divisions of genome at which to calculate ratios
	# Input 3: Length of genome
	# Output: Float of the maximum kernel density value of the 'hypothetical snps' ratio vector
	def self.max_hyp(hm, ht, div, genome_length)
		ratios = FitnessScore.ratio(hm, ht, div, genome_length)
		hyp = SNPdist.hyp_snps(ratios, genome_length)
		return max_density(hyp)
	end

	# Input 0: Array of homozygous snp positions
	# Input 1: Array of heterozygous snp positions
	# Input 2: Number of divisions of genome at which to calculate ratios
	# Input 3: Length of genome
	# Output: Integer of the total distance in bases, between adjacent SNPs from hypothetical ratio distribution (see SNPdist.hyp_snps)
	def self.hyp_distance(hm, ht, div, genome_length)
		ratios = FitnessScore.ratio(hm, ht, div, genome_length)
		hyp = SNPdist.hyp_snps(ratios, genome_length)
		return snp_distance(hyp)
	end
end