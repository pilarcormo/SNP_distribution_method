#encoding: utf-8
require_relative '../lib/model_genome'
require 'test/unit'

class TestModelGenome < Test::Unit::TestCase

	def test_fasta_to_char_array
		assert_equal(%w(A T G C A T A A A A A), ModelGenome::fasta_to_char_array('test/test/dummy.fasta'))
	end

	def test_get_snps
		hm_r = 'hm <- rnorm(150, 3000, 300)' 
		ht_r = 'ht <- runif(150, 1, 6000)'   
		hm, ht = ModelGenome::get_snps(hm_r, ht_r)
		assert(150 >= hm.length)
		assert(150 >= ht.length)
		assert_equal(hm.uniq.length, hm.length, 'hm not unique')
		assert_equal(ht.uniq.length, ht.length, 'ht not unique')
	end

	def test_get_frags
		seq = ['A']*1000
		frags = ModelGenome::get_frags(seq, 26)
		assert_equal(seq.length, frags.flatten.length, 'Extra frags!')
		assert_equal(seq, frags.flatten, 'Wrong bases!')
	end

	def test_pos_each_frag
		snp_pos = [1,4,6,15,18] # genome length = 20bp
		seq = ['A']*20
		frags = ModelGenome::get_frags(seq, 2)
		all_snps = ModelGenome::pos_each_frag(snp_pos, frags)
		assert_equal(snp_pos, all_snps[1].flatten, 'SNP positions in sub arrays not same')
		assert_equal(snp_pos.length, all_snps[0].flatten.length, 'Number of SNPs no right for positions on frags')
	end

	def test_fasta_array
		frags = [['A','T'], ['C','G'], ['A','T','C','G']]
		fasta_n_ids = ModelGenome::fasta_array(frags)
		assert_equal(['>frag1 Length = 2', 'AT'], fasta_n_ids[0])
	end

	def test_vcf_array
		frags = [['A','T'], ['C','G'], ['A','T','C','G']]
		pos_on_frags = [[1],[0],[1,2]]
		snp_pos_all = [1,2,5,6]
		hm = [2,5]
		ht = [1,6]
		vcf_array = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)
		assert_equal('frag1	1	.	T	A	100	PASS	AF=0.5', vcf_array[3])
		assert_equal('frag2	0	.	C	G	100	PASS	AF=1.0', vcf_array[4])
		assert_equal('frag3	2	.	C	G	100	PASS	AF=0.5', vcf_array[6])
	end
end