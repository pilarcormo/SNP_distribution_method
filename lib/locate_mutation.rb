#encoding: utf-8

class LocateMutation
	require 'rinruby'

	# Input 0: List of SNP positions
	# Input 1: The number of equally spaced points at which the density is to be estimated. Specify n as a power of two.
	# Output: The highest kernel density value for this SNP distribution
	def self.find_peak(snps, n)
		myr = RinRuby.new(echo=false)
		myr.n = n
		myr.snps = snps
		myr.eval 'kernel_density <- density(snps, n=n)'
		myr.eval 'index <- match(max(kernel_density$y),kernel_density$y)' # this only finds the first index with the max density if there is > 1
		myr.eval 'peak <- kernel_density$x[index]'
		peak = myr.pull 'peak'
		myr.quit
		return peak
	end

	

	# Input 0: Value under distribution peak (genome position as a float)
	# Input 1: List of homozygous SNP positions
	# Output: The closest homozygous SNP to the peak
	def self.closest_snp(peak, hm)
		hm.min{|a,b|  (peak-a).abs <=> (peak-b).abs }
	end
end