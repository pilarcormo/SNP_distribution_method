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



	def test_positions_by_fragment
		dic = {"frag1" => 1.0, "frag2"=>1.0, "frag3"=>2.0, "frag4"=>0.0}
		snp_list = [15, 18, 20, 25]
		assert_kind_of(Hash, dic)
		assert_kind_of(Array, snp_list)
		dic = Stuff.positions_by_fragment(dic, snp_list)
		assert_equal(dic, {"frag1" =>[15], "frag2"=>[18], "frag3"=>[20, 25]})
	end 

end 

