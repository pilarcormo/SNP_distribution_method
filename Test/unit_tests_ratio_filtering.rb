require 'test/unit'
require_relative '../lib/ratio_filtering'
require_relative '../lib/stuff'

class TestRatio_filtering < Test::Unit::TestCase
	def test_important_ratios 
		snp_hm = [0, 14, 20, 2]
		snp_ht = [5, 4, 2, 5]
		threshold = 0 
		adjust = 1
		ids = ["frag1", "frag2", "frag3", "frag4"]
		dic_ratios, ratios = Ratio_filtering.important_ratios(snp_hm, snp_ht, ids, threshold, adjust)
		assert_kind_of(Hash, dic_ratios)
		assert_kind_of(Array, ratios)
		assert_equal(dic_ratios, {"frag1" => 1/6.to_f, "frag2" => 3.to_f, "frag3" => 7.to_f, "frag4" => 1/2.to_f})
		assert_equal(ratios, [1/6.to_f, 3.to_f, 7.to_f, 1/2.to_f])
	end 
	def test_important_ids
		list = ["1", "3", "2", "4", "6", "7"]
		long_list = ["1", "5", "7", "36", "2", "3", "35", "20", "17"]
		short_list = Ratio_filtering.important_ids(list, long_list)
		assert_kind_of(Array, short_list)
		assert_equal(["1", "7", "2", "3"] ,short_list)
	end 
	def test_important_pos
		list = ["frag1", "frag2"]
		pos = {"frag1"=>[12, 13, 14], "frag2"=>[25], "frag4"=>[45]}
		pos_short = Ratio_filtering.important_pos(list, pos)
		assert_kind_of(Array, pos_short)
		assert_equal([12, 13, 14, 25], pos_short)
	end 