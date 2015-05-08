#encoding: utf-8

class SNPdist
	### Hypothetical SNP positions ###
	# Input 0: The contig size
	# Input 1: Array of ratios for homozygous to heterozygous SNPs at each contig. 
	# Output: A list of positions for the given ratios array size. 
	def self.general_positions(average_contig_size, ratios)
		positions = []
		ratios.length.times do |i|
			positions << average_contig_size
		end 
		x = 1 
		ratios.length.times do |i|
			positions[x] = positions[x-1].to_i + positions[x].to_i 
			x += 1
		end 
		return positions
	end 
	# Input 0: The homozygous, heterozygous or ratio densities in the ordered genome. 
	# Input 1: Previous output. A list of "hypothetical"  positions 
	# Output: A list of "hypothetical" SNP positions which represents the distribution of hm, ht SNPs or  homozygous/heterozygous SNP density ratio
	def self.densities_pos(raw_densities, positions)
		y = 0 
		densities_in_pos = []
		raw_densities.each do |density|
			(density*10).to_i.times do |i|
				densities_in_pos << positions[y]
			end 
			y += 1
		end 
		return densities_in_pos
	end 

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


	def self.plot_snps(snp_pos, correct_snps, location, genome_length, title, ylim)
		myr = RinRuby.new(echo = false)
		myr.snp_pos = snp_pos
		myr.correct_snps = correct_snps
		myr.genome_length = genome_length
		myr.location = location
		myr.title = title
		myr.ylim = ylim
		myr.eval 'png(paste(location,"/", "ratios",".png", sep=""))
		plot((1:512)*(genome_length/512), density(snp_pos)$y, xlim=c(0,genome_length), ylim=c(0,ylim), xlab="Genome length",
			ylab="Density", main=title)
		lines((1:512)*(genome_length/512), density(correct_snps)$y, lwd=3, col="#000099")
		dev.off()'
		myr.quit
	end


end