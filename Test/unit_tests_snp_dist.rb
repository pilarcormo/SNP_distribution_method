#encoding: utf-8
require_relative '../lib/snp_dist'
require_relative '../lib/write_it'
require 'test/unit'

class TestSNPdist < Test::Unit::TestCase

	def setup
		@ratios = WriteIt.file_to_floats_array("test/ratios_example.txt")
		@contig_size = 10
	end

	def test_general_positions
		pos = SNPdist.general_positions(@contig_size, @ratios)
		assert_equal([10, 20, 30, 40, 50, 60, 70, 70], pos)
	end

	def test_densities_pos
		@pos = [10, 20, 30, 40, 50, 60, 70]
		densities = SNPdist.densities_pos(@ratios, @pos)
		assert_equal([10, 20, 30, 30, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 70], densities)
	end 

	def test_ylim
		ylim = SNPdist.get_ylim([1,2,3,4,5], 2000, 'density')
		assert_equal(0.25, ylim)
	end
end