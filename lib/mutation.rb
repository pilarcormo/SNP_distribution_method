#encoding: utf-8
require_relative 'locate_mutation'
require_relative 'snp_dist'
require_relative 'reform_ratio'
require_relative 'plot'
require_relative 'write_it'
require 'pp'

class Mutation
	def self.test_genomes_define(hm, ht, perm_hm, perm_ht, genome_length, ratios, expected_ratios) 
		hm.pop
		div = 100
		n = 2000
		hyp = SNPdist.hyp_snps(ratios, genome_length)
		peak =  LocateMutation.find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
		causal = LocateMutation.closest_snp(peak, hm)
		perm_hyp = SNPdist.hyp_snps(expected_ratios, genome_length)

		perm_peak = LocateMutation.find_peak(perm_hyp, n)
		candidate = LocateMutation.closest_snp(perm_peak, perm_hm)
		normalised = (candidate - causal).abs
		percent = (normalised*100)/genome_length.to_f
		return causal, candidate, percent
	end

	def self.test_genomes_distribution_plot(genome_length, ratios, expected_ratios, dataset, perm)
		hm, ht, hyp, ylim_hm, ylim_ht, ylim_hyp = [],[],[],[],[],[]
		Dir.chdir(File.join(Dir.home, "SNP_distribution_method/Small_genomes/#{perm}")) do
			hom_snps = WriteIt.file_to_ints_array("hm_snps_short.txt")
			hm << hom_snps
			ylim_hm << SNPdist.get_ylim(hom_snps, genome_length, 'density')

			het_snps = WriteIt.file_to_ints_array("ht_snps_short.txt")
			ht << het_snps
			ylim_ht << SNPdist.get_ylim(het_snps, genome_length, 'density')

			hyp_snps = SNPdist.hyp_snps(expected_ratios, genome_length)
			hyp << hyp_snps
			ylim_hyp << SNPdist.get_ylim(hyp_snps, genome_length, 'density')
		end

		Dir.chdir(File.join(Dir.home, "SNP_distribution_method/Small_genomes/#{perm}")) do

			perm_hm = WriteIt.file_to_ints_array("perm_hm.txt")
			SNPdist.plot_snps(perm_hm, hm[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'hm',
				'Homozygous SNP density', ylim_hm[0])

			perm_ht = WriteIt.file_to_ints_array("perm_ht.txt")
			SNPdist.plot_snps(perm_ht, ht[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'ht',
				'Heterozygous SNP density', ylim_ht[0])

			perm_hyp = SNPdist.hyp_snps(ratios, genome_length)
			SNPdist.plot_snps(perm_hyp, hyp[0], "SNP_distribution_method/Small_genomes", "#{perm}", 1, genome_length, 'hyp', 
				'Approximated ratio of homozygous to heterozygous SNP density', ylim_hyp[0])
		end
	end


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