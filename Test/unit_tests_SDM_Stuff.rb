#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/SDM'
require_relative '../lib/reform_ratio'
require_relative '../lib/stuff'

#Take 4 numbers generated randomly that follow a normal distribution, rnorm(4, mean=10, sd=2)

# hash = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}

class TestSDM < Test::Unit::TestCase
	def setup 
		@dic_hm_inv = {11.374979=>["a"], 11.842904=>["b"], 10.120768=>["c"], 5.022447=>["d", "e", "f"]}
	end
	def test_sorting
		perm = SDM::sorting(@dic_hm_inv)
		assert_equal(["e", "d", "a", "b", "c", "f"], perm)
	end
end 

class TestSuff < Test::Unit::TestCase
	def setup
		@vcf_file = "test/test.vcf"
		@fasta_file = "test/test.fasta"
		@f_array = ReformRatio::fasta_array("test/test.fasta")
		# @snp_data = ReformRatio::get_snp_data('arabidopsis_datasets/dataset_small2kb/snps.vcf')
		# @genome_length = ReformRatio::genome_length(@fasta_array)
	end
	def test_snps_in_vcf
		hm, ht = Stuff.snps_in_vcf(@vcf_file)
		assert_equal(["frag1", "frag1"], hm)
		assert_equal(["frag2", "frag3"], ht)
	end 
	def test_create_hash_snps
		hm = ["frag1", "frag1"]
		ht = ["frag2", "frag3"]
		dic_hm, dic_ht = Stuff.create_hash_snps(hm, ht)
		assert_kind_of(Hash, dic_hm)
		assert_kind_of(Hash, dic_ht)
		assert_equal(dic_hm, {"frag1"=>"2"})
		assert_equal(dic_ht, {"frag2"=>"1", "frag3"=>"1"})
	end
	def test_define_snps
		ids = ["frag1", "frag3", "frag2"]
		dic1 = {"frag1"=>"2"}
		dic2 = {"frag2"=>"1", "frag3"=>"1"}
		shu_dic1, shu_dic2, snps_1, snps_2 = Stuff.define_snps(ids, dic1, dic2)
		assert_kind_of(Hash, shu_dic1)
		assert_kind_of(Hash, shu_dic2)
		assert_kind_of(Array, snps_1)
		assert_kind_of(Array, snps_2)
		assert_kind_of(Float, snps_1[0])
		assert_equal(shu_dic1, {"frag1"=>2.0, "frag3"=>0, "frag2"=>0})
		assert_equal(shu_dic2, {"frag1"=>0, "frag3"=>1.0, "frag2"=>1.0})
		assert_equal(snps_1, [2.0, 0.0, 0.0])
		assert_equal(snps_2, [0, 1.0, 1.0])
	end 
	def test_normalise_by_length
		dic1 = {"frag1"=>2.0, "frag3"=>0, "frag2"=>0}
		dic2 = {"frag1"=>0, "frag3"=>1.0, "frag2"=>1.0}
		snps_1 = [2.0, 0.0, 0.0]
		snps_2 = [0, 1.0, 1.0]
		ids, lengths = ReformRatio.fasta_id_n_lengths(@f_array)
		dic_norm1, dic_norm2 = Stuff.normalise_by_length(ids, dic1, dic2, lengths)
		assert_equal(dic_norm1, {"frag1"=>2.0/11, "frag3"=>0/8, "frag2"=>0/8})
		assert_equal(dic_norm2, {"frag1"=>0/11, "frag3"=>1.0/8, "frag2"=>1.0/8})

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
		ids = ["frag1", "frag2", "frag3", "frag4"]
		dic_ratios, ratios = Stuff.important_ratios(snp_hm, snp_ht, ids)
		assert_kind_of(Hash, dic_ratios)
		assert_kind_of(Array, ratios)
		assert_equal(dic_ratios, {"frag2" => 3.to_f, "frag3" => 7.to_f})
		assert_equal(ratios, [3.to_f, 7.to_f])
	end 
end 

