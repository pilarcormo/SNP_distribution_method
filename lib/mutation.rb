#encoding: utf-8
require_relative 'locate_mutation'
require_relative 'snp_dist'
require_relative 'reform_ratio'
require_relative 'plot'
require_relative 'write_it'
require 'pp'

class Mutation
	def self.candidate(mut, frag_pos_hm)
		candidate_mutations ={}
		mut.each do |frag|
			if frag_pos_hm.has_key?(frag)
				candidate_mutations.store(frag, frag_pos_hm[frag])
			end
		end 
		return candidate_mutations
	end 

	def self.density_plots(contig_size, ratios, expected_ratios, snps_hm, snps_ht, center, file, mut, frag_pos_hm) 
		n = 1048576*4
		average_positions = SNPdist.general_positions(contig_size, ratios)
		hyp_ratios = SNPdist.densities_pos(expected_ratios, average_positions)
		real_ratios = SNPdist.densities_pos(ratios, average_positions)
		WriteIt::write_txt("#{file}/hyp_ratios", hyp_ratios)
		WriteIt::write_txt("#{file}/ratios", real_ratios)
		peak =  LocateMutation.find_peak(hyp_ratios, n) # Find the peak in the approximated (hypothetical SNP) distribution
		ylim = Plot.get_ylim(hyp_ratios, center)
		candidate_peak = LocateMutation.closest_snp(peak, snps_hm)
		Plot::densities(snps_hm, snps_ht, hyp_ratios, center, file)
		Plot::comparison(real_ratios, hyp_ratios, center, file, ylim)
		Plot::qqplot(snps_hm, file, "QQplot for hm density", "Theoretical normal distribution", "Hypothetical SNP density")
		#Plot::qqplot(hyp_ratios, file, "QQplot for the ratios", "Theoretical normal distribution", "Hypothetical ratios")
		candidate_mutations = Mutation.candidate(mut, frag_pos_hm)
		File.open("#{file}/mutation.txt", "w+") do |f|
			f.puts "The length of the group of contigs that form the peak of the distribution is #{center.to_i} bp"
			f.puts "The mutation is likely to be found on the following contigs #{candidate_mutations}"
		end 
		return candidate_peak
	end
end