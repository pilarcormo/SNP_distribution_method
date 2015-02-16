#encoding: utf-8
require_relative '../lib/locate_mutation'
require 'test/unit'

class TestLocateMutation < Test::Unit::TestCase

	def setup
		@example_snps = [105,109,87,96,111,95,100,88,110,92]
	end

	def test_find_peak
		assert_equal(94, (LocateMutation.find_peak(@example_snps, 512)))
	end

	def test_closest
		assert_equal(95, LocateMutation.closest_snp(95.345, @example_snps))
		assert_equal(100, LocateMutation.closest_snp(101, @example_snps))
		assert_equal(111, LocateMutation.closest_snp(200, @example_snps))
	end
end