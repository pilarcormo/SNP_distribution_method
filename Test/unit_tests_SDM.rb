#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/SDM'
class TestSDM < Test::Unit::TestCase
	def setup 
		@dic_ratio_inv_back = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}
		@dic_hm_inv_back = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}
		@dic_hm_inv_out = {11.374979=>["a1", "a2", "a3", "a4", "a5"], 11.842904=>["b1", "b2", "b3", "b4"], 10.120768=>["c1", "c2"], 5.022447=>["d1", "e1", "f1", "d2", "e2", "f2", "d3", "e3", "f3", "f4"]}
		@dic_ratio_inv_out = {11.374979=>["a1", "a2", "a3", "a4", "a5"], 11.842904=>["b1", "b2", "b3", "b4"], 10.120768=>["c1", "c2"], 5.022447=>["d1", "e1", "f1", "d2", "e2", "f2", "d3", "e3", "f3", "f4"]}
		@cross_back = "back"
		@cross_out = "out"
		@dic_pos_hm_back = {"a"=>[1], "b"=>[2], "c"=>[3, 4], "d"=>[5, 6], "e"=>[7], "f"=>[10]}
		@dic_pos_hm_out = {"a1"=>[1], "a2"=>[2], "a3"=>[5], "a4"=>[6], "a5"=>[7], "b1"=>[8], "b2"=>[9], "b3"=>[10], "b4"=>[17], "c1"=>[21], "c2"=>[23], "d1"=>[25], "e1"=>[34], "f1"=>[37], "d2"=>[50], "e2"=>[51], "f2"=>[53], "d3"=>[56], "e3"=>[57], "f3"=>[60], "f4"=>[70]}

	end
	def test_sorting 
		right, left = [], []
		perm_back, mut_back = SDM::sorting(@dic_hm_inv_back, @cross_back)
		perm_out, mut_out = SDM::sorting(@dic_hm_inv_out, @cross_out)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_back)
		assert_equal([ "a", "b"], mut_back)
		assert_equal(["d1","e1", "f1", "d2", "e2", "c1", "a2", "a3", "a1", "b1", "b2", "b4", "b3", "a5", "a4", "c2", "f4", "f3", "e3", "d3", "f2"], perm_out)
		assert_equal(["e1", "f1", "d2", "e2", "c1", "a2", "a3", "a1", "b1", "b2", "b4", "b3", "a5", "a4", "c2", "f4", "f3", "e3", "d3", "f2"], mut_out)	
	end 

	def test_calling_SDM
		perm_hm, perm_ratio, mut, hyp = SDM::calling_SDM(@dic_hm_inv_back, @dic_ratio_inv_back, @cross_back, @dic_pos_hm_back)
		perm_hm_out, perm_ratio_out, mut_out, hyp_out = SDM::calling_SDM(@dic_hm_inv_out, @dic_ratio_inv_out, @cross_out, @dic_pos_hm_out)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_hm)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm_ratio)
		assert_equal([ "a", "b"], mut)
		assert_equal([1, 2], hyp)
		assert_equal(["d1","e1", "f1", "d2", "e2", "c1", "a2", "a3", "a1", "b1", "b2", "b4", "b3", "a5", "a4", "c2", "f4", "f3", "e3", "d3", "f2"], perm_hm_out)
		assert_equal(["d1","e1", "f1", "d2", "e2", "c1", "a2", "a3", "a1", "b1", "b2", "b4", "b3", "a5", "a4", "c2", "f4", "f3", "e3", "d3", "f2"], perm_ratio_out)
		assert_equal(["e1", "f1", "d2", "e2", "c1", "a2", "a3", "a1", "b1", "b2", "b4", "b3", "a5", "a4", "c2", "f4", "f3", "e3", "d3", "f2"], mut_out)
		assert_equal([1, 2, 5, 6, 7, 8, 9, 10, 17, 21, 23, 34, 37, 50, 51, 53, 56, 57, 60, 70], hyp_out)
	end 

end 
