#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/reform_ratio'
require_relative '../lib/stuff'

class TestSuff < Test::Unit::TestCase
	def setup
		@vcf_file = "test/test.vcf"
		@fasta_file = "test/test.fasta"
		@f_array = ReformRatio::fasta_array("test/test.fasta")
	end
	def test_snps_in_vcf
		snp_data, hm, ht = Stuff.snps_in_vcf(@vcf_file)
		assert_equal(["frag1", "frag1"], hm)
		assert_equal(["frag2", "frag3"], ht)
		assert_equal([["frag1", "frag1", "frag2", "frag3"], [7, 8, 2, 2], {"frag1" =>2, "frag2" =>1, "frag3" =>1}, [{"AF"=>"1.0"},{"AF"=>"1.0"},{"AF"=>"0.5"},{"AF"=>"0.5"}]], snp_data)
	end 
	def test_safe_invert
		hash = {"frag1" => 1, "frag2" => 1, "frag3" => 1, "frag4" => 2, "frag5" => 1}
		hash_inv = Stuff.safe_invert(hash)
		assert_equal({1 => ["frag1",  "frag2", "frag3", "frag5"], 2 => ["frag4"]}, hash_inv)
	end 
	def test_dic_id_pos
		hm = ["frag1", "frag1", "frag1", "frag4"]
		ht = ["frag1", "frag1", "frag2", "frag4"]
		pos1 = [12, 13, 14, 45]
		pos2 = [1, 4, 25, 40]
		dic_pos_hm = Stuff.dic_id_pos(hm, pos1)
		dic_pos_ht = Stuff.dic_id_pos(ht, pos2)
		assert_kind_of(Hash, dic_pos_hm)
		assert_kind_of(Hash, dic_pos_ht)
		assert_equal(dic_pos_hm, {"frag1"=>[12, 13, 14], "frag4"=>[45]})
		assert_equal(dic_pos_ht, {"frag1"=> [1, 4], "frag2"=>[25], "frag4"=>[40]})
	end
	def test_create_hash_number
		hm = ["frag1", "frag1"]
		ht = ["frag2", "frag3"]
		dic_hm = Stuff.create_hash_number(hm)
		dic_ht = Stuff.create_hash_number(ht)
		assert_kind_of(Hash, dic_hm)
		assert_kind_of(Hash, dic_ht)
		assert_equal(dic_hm, {"frag1"=>"2"})
		assert_equal(dic_ht, {"frag2"=>"1", "frag3"=>"1"})
	end
	def test_define_snps
		ids = ["frag1", "frag3", "frag2"]
		dic1 = {"frag1"=>"2"}
		shu_dic1, snps_1  = Stuff.define_snps(ids, dic1)
		assert_kind_of(Hash, shu_dic1)
		assert_kind_of(Array, snps_1)
		assert_kind_of(Float, snps_1[0])
		assert_equal(shu_dic1, {"frag1"=>2.0, "frag3"=>0, "frag2"=>0})
		assert_equal(snps_1, [2.0, 0.0, 0.0])
	end 
	def test_normalise_by_length ##problem with float, check original 
		dic1 = {"frag1"=>2.0, "frag3"=>0, "frag2"=>0}
		snps_1 = [2.0, 0.0, 0.0]
		ids, lengths = ReformRatio.fasta_id_n_lengths(@f_array)
		dic_norm1 = Stuff.normalise_by_length(lengths, dic1)
		assert_equal(dic_norm1, {2.0/11.to_f =>["frag1"], 0/8.to_f =>["frag3","frag2"]})

	end
	def test_create_perm_fasta
		perm = []
		@fasta_array = ReformRatio::fasta_array("test/test2.fasta")
		ids, lengths = ReformRatio.fasta_id_n_lengths(@fasta_array)
		perm = ["frag1", "frag3", "frag2"]
		fasta_perm = Stuff.create_perm_fasta(perm, @fasta_array, ids)
		assert_equal(fasta_perm, [">frag1 Length = 8", "CCAAATAC\n", ">frag3 Length = 7", "ACGACAC\n", ">frag2 Length = 8", "GCAATCGG\n"])
	end 

	def test_positions_by_fragment
		dic = {"frag1" => 1.0, "frag2"=>1.0, "frag3"=>2.0, "frag4"=>0.0}
		snp_list = [15, 18, 20, 25]
		assert_kind_of(Hash, dic)
		assert_kind_of(Array, snp_list)
		dic = Stuff.positions_by_fragment(dic, snp_list)
		assert_equal(dic, {"frag1" =>[15], "frag2"=>[18], "frag3"=>[20, 25]})
	end 
	def test_important_ratios 
		snp_hm = [0, 14, 20, 2]
		snp_ht = [5, 4, 2, 5]
		threshold = 0 
		adjust = 1
		ids = ["frag1", "frag2", "frag3", "frag4"]
		dic_ratios, ratios = Stuff.important_ratios(snp_hm, snp_ht, ids, threshold, adjust)
		assert_kind_of(Hash, dic_ratios)
		assert_kind_of(Array, ratios)
		assert_equal(dic_ratios, {"frag1" => 1/6.to_f, "frag2" => 3.to_f, "frag3" => 7.to_f, "frag4" => 1/2.to_f})
		assert_equal(ratios, [1/6.to_f, 3.to_f, 7.to_f, 1/2.to_f])
	end 
	def test_important_ids
		list = ["1", "3", "2", "4", "6", "7"]
		long_list = ["1", "5", "7", "36", "2", "3", "35", "20", "17"]
		short_list = Stuff.important_ids(list, long_list)
		assert_kind_of(Array, short_list)
		assert_equal(["1", "7", "2", "3"] ,short_list)
	end 
	def test_important_pos
		list = ["frag1", "frag2"]
		pos = {"frag1"=>[12, 13, 14], "frag2"=>[25], "frag4"=>[45]}
		pos_short = Stuff.important_pos(list, pos)
		assert_kind_of(Array, pos_short)
		assert_equal([12, 13, 14, 25], pos_short)
	end 
end 

