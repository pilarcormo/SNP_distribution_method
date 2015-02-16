#encoding: utf-8
# require_relative 'lib/write_it'
require_relative 'locate_mutation'
require_relative 'snp_dist'
require_relative 'reform_ratio'
require_relative 'fitness_score'
require 'csv'

class Mutation
	def self.define(hm, ht, perm_hm, perm_ht, genome_length, ratios, expected_ratios) 
		div = 100
		n = 1048576*4
		# hyp = SNPdist.hyp_snps(ratios, genome_length)
		# puts hyp.length
		

		peak =  LocateMutation.find_peak(hm, n) # Find the peak in the approximated (hypothetical SNP) distribution
 
		causal = LocateMutation.closest_snp(peak, hm)

		# Dir.chdir(File.join(Dir.home, "small_genomes_SNPs/arabidopsis_datasets")) do
		# 	# CSV.open("ratios.csv", "ab") do |csv|
		# 	# 	csv << ["ratios", "perm_ratios"]
		# 	# end
		# 	ratios.each do |i|
		# 		CSV.open("expected_ratios.csv", "ab") do |csv|
		# 			csv << [i]
		# 		end
		# 	end
		# 	expected_ratios.each do |o|
		# 		CSV.open("expected_ratios.csv", "ab") do |csv|
		# 			csv << [o]
		# 		end
		# 	end 
		# end

		# perm_hyp = SNPdist.hyp_snps(expected_ratios, genome_length)


		perm_peak = LocateMutation.find_peak(perm_hm, n)
		candidate = LocateMutation.closest_snp(perm_peak, perm_hm)

		puts "Location of causal mutation in correctly ordered genome: #{causal}"
		puts "Candidate SNP position in permutation: #{candidate}"

		normalised = (candidate - causal).abs
		
		percent = (normalised*100)/genome_length.to_f
		puts "Shift #{percent} %"

		return percent
	end
end