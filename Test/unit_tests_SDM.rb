#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/SDM'
class TestSDM < Test::Unit::TestCase
	def setup 
		@f_array = ReformRatio::fasta_array("test/test.fasta")
		@dic_ratio_inv_back = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}
		@dic_hm_inv_back = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}
		@dic_pos_hm_back = {"a"=>[1], "b"=>[2], "c"=>[3, 4], "d"=>[5, 6], "e"=>[7], "f"=>[10]}
	end
	def test_normalise_by_length 
		dic1 = {"frag1"=>2.0, "frag3"=>0, "frag2"=>0}
		snps_1 = [2.0, 0.0, 0.0]
		ids, lengths = ReformRatio.fasta_id_n_lengths(@f_array)
		dic_norm1 = SDM.normalise_by_length(lengths, dic1)
		assert_equal(dic_norm1, {2.0/11.to_f =>["frag1"], 0/8.to_f =>["frag3","frag2"]})
	end

	def test_sorting 
		right, left = [], []
		perm_back, mut_back = SDM::sorting(@dic_hm_inv_back)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_back)
		assert_equal(["e", "d", "a", "b", "c", "f"], mut_back)
	end 

	def test_calling_SDM
		perm_hm, perm_ratio, mut, hyp = SDM::calling_SDM(@dic_hm_inv_back, @dic_ratio_inv_back, @dic_pos_hm_back)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_hm)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_ratio)
		assert_equal(["e", "d", "a", "b", "c", "f"], mut)
		assert_equal([1, 2, 3, 4, 5, 6, 7, 10], hyp)
	end 

end 
