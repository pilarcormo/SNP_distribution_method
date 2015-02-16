#encoding: utf-8
require_relative '../lib/reform_ratio'
require 'test/unit'

class FakeFasta
	attr_accessor :entry_id, :length
	def initialize
		@entry_id = 'FragX'
		@length = 10
	end
end

class TestReform < Test::Unit::TestCase

	def setup
		vcf_file = "test/test/dummy.vcf"
		fasta_file = "test/test/dummy.fasta"

		snp_data = ReformRatio::get_snp_data(vcf_file)
		fasta = ReformRatio::fasta_array(fasta_file)

		snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) # array of no. of snps per frag in same order as fasta
		@pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) # get snp positions for each frag in array of arrays
		@actual_pos = ReformRatio::total_pos(@pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
		@het_snps, @hom_snps = ReformRatio::het_hom(@actual_pos, @pos_n_info[1])

		@perm_pos = ReformRatio::perm_pos(fasta, snp_data)
	end

	def fasta
		frag1, frag2, frag3 = FakeFasta.new, FakeFasta.new, FakeFasta.new
		frag1.entry_id, frag2.entry_id, frag3.entry_id = 'frag1', 'frag2', 'frag3'
		[frag1, frag2, frag3]
	end

	def test_total_pos
		pos = [[1,2,3],[4,5,6],[7,8,9]]
		lengths = [10,10,10]
		expected = [1,2,3,14,15,16,27,28,29]
		assert_equal(expected, ReformRatio::total_pos(pos, lengths))

		assert_equal([7,8,13,21], @actual_pos)
	end

	def test_fasta_id_n_lengths
		ids_n_lengths = ReformRatio::fasta_id_n_lengths(fasta)
		ids = ids_n_lengths[0]
		lengths = ids_n_lengths[1]
		assert_equal(['frag1', 'frag2', 'frag3'], ids)
		assert_equal([10,10,10], lengths)
	end

	def test_snps_per_fasta_frag
		h = {'frag2'=>3,'frag1'=>2,'frag3'=>5}
		assert_equal([2,3,5], ReformRatio::snps_per_fasta_frag(h, fasta))
	end

	def test_get_positions
		vcfs_chrom = ['frag3', 'frag3', 'frag3', 'frag1', 'frag1', 'frag2']
		vcfs_pos = [2,3,5,4,7,5] # f3= 2,3,5  f1= 4,7  f2= 5
		vcfs_info = [{'AF'=>'snp1'},{'AF'=>'snp2'},{'AF'=>'snp3'},{'AF'=>'snp4'},{'AF'=>'snp5'},{'AF'=>'snp6'}] # f3= s1,2,3  f1= s4,5  f2= s6
		snps_per_frag = [2,1,3] # in same order as fasta array
		pos_n_info = ReformRatio::get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag, vcfs_info)
		pos = pos_n_info[0]
		info = pos_n_info[1]
		assert_equal([[4,7],[5],[2,3,5]], pos)
		assert_equal([[{'AF'=>'snp4'},{'AF'=>'snp5'}] ,[{'AF'=>'snp6'}], [{'AF'=>'snp1'},{'AF'=>'snp2'},{'AF'=>'snp3'}]], info)

		assert_equal([[7,8],[2],[2]], @pos_n_info[0])
		assert_equal([[{'AF'=>'1.0','NS'=>'5'},{'AF'=>'1.0'}], [{'AF'=>'0.5'}], [{'AF'=>'1.0'}]], @pos_n_info[1])
	end

	def test_fasta_array
		fasta_array = ReformRatio::fasta_array('test/test/dummy.fasta')
		assert_equal('frag1', fasta_array[0].entry_id)
		assert_equal('AAAAAAAA', fasta_array[1].seq)
		assert_equal(8, fasta_array[2].length)
	end

	def test_get_snp_data
		vcfs_chrom = %w(frag1 frag1 frag2 frag3)
		vcfs_pos = [7,8,2,2]
		num_snps_frag_hash = {'frag1'=>2, 'frag2'=>1, 'frag3'=>1}
		vcfs_info = [{'AF'=>'1.0','NS'=>'5'}, {'AF'=>'1.0'}, {'AF'=>'0.5'}, {'AF'=>'1.0'}]
		snp_data = [vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info]
		assert_equal(snp_data, ReformRatio::get_snp_data('test/test/dummy.vcf'))
	end

	def test_het_hom
		vcfs_info = [{'AF'=>'1.0'}, {'AF'=>'1.0'}, {'AF'=>'0.5'}, {'AF'=>'0.5'}, {'AF'=>'1.0'}]
		actual_pos = [2, 17, 56, 190, 191]
		hom = [2, 17, 191]
		het = [56, 190]
		assert_equal([het,hom], ReformRatio::het_hom(actual_pos, vcfs_info))

		assert_equal([13], @het_snps)
		assert_equal([7,8,21], @hom_snps)
	end

	def test_perm_pos
		ht, hm = @perm_pos
		assert_equal([13], ht)
		assert_equal([7,8,21], hm)
	end
end

