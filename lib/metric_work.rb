#encoding: utf-8

class MetricWork

	require_relative 'reform_ratio'
	require_relative 'GATOC'
	require 'pmeth'
	require 'pdist'

=begin Testing how well the fitness method of the genetic algorithm,
		identifies permuations that are approaching the correct order.

		1. Get an array of correctly ordered FASTA contigs
		2. Create permutations that are progressively further from the correct
		3. Save the permutation fitness scores directly into a csv, to make a ggplot (but can also save permutation txt files - possibly not neccesary)
		4. Create a plot to represent this
=end

	def self.fitness(fasta, snp_data, genome_length, method, params)
		opts = {:div => nil, :expected_ratios => nil}.merge!(params)
		het_snps, hom_snps = ReformRatio.perm_pos(fasta, snp_data)
		case method
		when 'count_ratio' then score = FitnessScore.count_ratio(hom_snps, het_snps, opts[:div], genome_length, opts[:expected_ratios])
		when 'snp_distance' then score = FitnessScore.snp_distance(hom_snps)
		when 'max_density' then score = FitnessScore.max_density(hom_snps)
		when 'max_ratio' then score = FitnessScore.max_ratio(hom_snps, het_snps)
		when 'max_hyp' then score = FitnessScore.max_hyp(hom_snps, het_snps, opts[:div], genome_length)
		when 'hyp_distance' then score = FitnessScore.hyp_distance(hom_snps, het_snps, opts[:div], genome_length)
		end
		return score.to_f
	end

	def self.score(fasta, perm, snp_data, genome_length, metric, params)
		case metric
		when 'DeviationDistance'
			score = PDist.deviation(fasta, perm)
		when 'SquareDeviationDistance'
			score = PDist.square(fasta, perm)
		when 'HammingDistance'
			score = PDist.hamming(fasta, perm)
		when 'RDistance'
			score = PDist.rdist(fasta, perm)
		when 'LongestCommonSubsequence'
			score = PDist.lcs(fasta, perm)
		when 'KendallsTau'
			score = PDist.kendalls_tau(fasta, perm)
		else
			score = fitness(perm, snp_data, genome_length, metric.gsub(/(?<=[a-z])(?=[A-Z])/, ' ').downcase.tr(' ', '_'), params) # fitness method
		end
		return score
	end

	def self.adjacent_swaps_csv(dataset, size, pop_num, metric, filename, swap_num, params)
		# 1
		fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta") # correct permutation
		###

		# 2
		snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{dataset}/snps.vcf")
		genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")
		###

		# 3/4: Population of adjacent_swap mutants (of the correct contig order)
		start_pop = []
		size.times do
			start_pop << fasta
		end

		shuf_scores = []
		size.times do
			shuf_scores << score(fasta, fasta.shuffle, snp_data, genome_length, metric, params)
		end

		WriteIt.add_to("arabidopsis_datasets/#{dataset}/#{filename}.csv", "population,#{metric},shuffled")
		x = y = 1
		pop_num.times do
			adj_pop = []
			start_pop.each do |perm|
				new_perm = perm
				swap_num.times do
					new_perm = PMeth.adjacent_swap(new_perm) 
				end
				puts "adjacent_swaps: another #{swap_num} for pop: #{x}"
				adj_pop << new_perm # need this population to be the next starting population
				score = MetricWork.score(fasta, perm, snp_data, genome_length, metric, params)
				if y <= size
					WriteIt.add_to("arabidopsis_datasets/#{dataset}/#{filename}.csv", "#{x*swap_num},#{score},#{shuf_scores[y-1]}")
				else
					WriteIt.add_to("arabidopsis_datasets/#{dataset}/#{filename}.csv", "#{x*swap_num},#{score}")
				end
				y+=1
			end
			start_pop = adj_pop
			x+=1
		end
	end

	# 5:
	def self.metric_test_plot(dataset, filename, metric, x_axis, y_axis, title, input_file)
		myr = RinRuby.new(echo = false)
		myr.dataset = dataset
		myr.filename = filename
		myr.metric = metric
		myr.title = title
		myr.x_axis = x_axis
		myr.y_axis = y_axis
		myr.input_file = input_file
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		myr.eval "df <- read.csv(paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset, '/', input_file, '.csv', sep=''))"
		myr.eval "p <- metric_test_plot(df, title, x_axis, y_axis, metric)"
		myr.eval "ggsave(p, file = paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset,'/', filename,'.png', sep = ''))"
		myr.quit
		puts 'made a plot'
	end

	
end
