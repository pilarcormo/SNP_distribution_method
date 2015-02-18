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

		hyp = SNPdist.hyp_snps(ratios, genome_length)
	
		peak =  LocateMutation.find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
 
		causal = LocateMutation.closest_snp(peak, hm)
		perm_hyp = SNPdist.hyp_snps(expected_ratios, genome_length)
		perm_peak = LocateMutation.find_peak(perm_hyp, n)
		candidate = LocateMutation.closest_snp(perm_peak, perm_hm)

		puts "Location of causal mutation in correctly ordered genome: #{causal}"
		puts "Candidate SNP position in permutation: #{candidate}"

		normalised = (candidate - causal).abs
		
		percent = (normalised*100)/genome_length.to_f
		puts "Shift #{percent} %"

		return percent
	end
	def self.distribution_plot(hom_snps, het_snps, perm_hm, perm_ht, genome_length, ratios, expected_ratios)
		hm, ht, hyp, ylim_hm, ylim_ht, ylim_hyp = [],[],[],[],[],[]
		hm << hom_snps
		ylim_hm << SNPdist.get_ylim(hom_snps, genome_length, 'density')
		ht << het_snps
		ylim_ht << SNPdist.get_ylim(het_snps, genome_length, 'density')
		hyp_snps = SNPdist.hyp_snps(expected_ratios, genome_length)
		hyp << hyp_snps
		ylim_hyp << SNPdist.get_ylim(hyp_snps, genome_length, 'density')
		plot_hm = SNPdist.plot_snps(perm_hm, hm[0], "SNP_distribution_method/arabidopsis_datasets", "dataset_small2kb/1802", genome_length, 'hm', 'Homozygous SNP density', ylim_hm[0])
		plot_ht = SNPdist.plot_snps(perm_ht, ht[0], "SNP_distribution_method/arabidopsis_datasets", "dataset_small2kb/1802", genome_length, 'ht', 'Heterozygous SNP density', ylim_ht[0])
		perm_hyp = SNPdist.hyp_snps(ratios, genome_length)
		plot_hyp = SNPdist.plot_snps(perm_hyp, hyp[0], "SNP_distribution_method/arabidopsis_datasets", "dataset_small2kb/1802", genome_length, 'hyp', 'Approximated ratio of homozygous to heterozygous SNP density', ylim_hyp[0])
		return plot_hm, plot_ht, plot_hyp
	end
end