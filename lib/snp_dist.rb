#encoding: utf-8

class SNPdist
	require 'rinruby'

	### Hypothetical SNP positions ###

	# Input 0: Array of ratios for homozygous to heterozygous SNPs at divisions of the genome
	# Input 1: The length of the genome
	# Output: A list of "hypothetical" SNP positions which represents the distribution of homozygous/heterozygous SNP density ratio
	def self.hyp_snps(ratios, genome_length)
		breaks = []
		(1..ratios.length).to_a.each do |i|
			breaks << (genome_length/ratios.length.to_f)*i
		end
		# pp "This are breaks #{breaks}"
		hyp, x = [], 0
		ratios.each do |ratio| 
			(ratio*10).to_i.times do
				hyp << rand(genome_length/ratios.length.to_f) + breaks[x] # random value from within the range that the freq has been taken
			end
			x+=1
		end
		# pp "hyp hyp #{hyp}"
		return hyp # These don't need to be unique or integers like the real SNPs, since they are just representing a distribution
	end

	### Plotting Methods ###

	# Input 0: An array that contains ratio values for homozgous to heterozygous at divisions of a genome
	# Input 1: The length of the genome
	# Input 2: A string of either 'density' or ratio
	# Output: Float of the highest value on the y axis when the ratio is plotted in the two different ways shown in the method
	def self.get_ylim(array, genome_length, plot_type)
		myr = RinRuby.new(echo = false)
		myr.assign 'array', array
		myr.assign 'genome_length', genome_length
		if plot_type == 'density'
			myr.eval 'plot((1:512)*(genome_length/512), density(array)$y)'
		end
		ylim = myr.pull 'par("yaxp")[2] + par("yaxp")[2]/par("yaxp")[3]'
		myr.quit
		return ylim
	end

	# Input 0: A list of SNP positions for the permutation
	# Input 1: A list of SNP positions for the correct permutation
	# Input 1: Location at which to save the plot
	# Input 2: The dataset to use (and the sub-folder to save the plot in)
	# Input 3: The generation of the genetic algorithm the ratio is being plotted for
	# Input 4: The length of the genome
	# Input 5: String indicating the type of SNPs that Input 0 are
	# Input 6: Title of plot
	# Input 7: The highest value on the y axis for the plot
	# Output: Plot of kernel density estimate for the SNPs over the genome
	def self.plot_snps(snp_pos, correct_snps, location, dataset_run, gen, genome_length, type, title, ylim)
		myr = RinRuby.new(echo = false)
		myr.snp_pos = snp_pos
		myr.location = location
		myr.correct_snps = correct_snps
		myr.genome_length = genome_length
		myr.dataset_run = dataset_run
		myr.gen = gen
		myr.type = type
		myr.title = title
		myr.ylim = ylim
		myr.eval 'png(paste("~/",location,"/", dataset_run, "/Plot_", type, ".png", sep=""))
		plot((1:512)*(genome_length/512), density(snp_pos)$y, xlim=c(0,genome_length), ylim=c(0,ylim), xlab="Genome length",
			ylab="Density", main=title)
		lines((1:512)*(genome_length/512), density(correct_snps)$y, lwd=3, col="#000099")
		dev.off()'
		myr.quit
	end
end