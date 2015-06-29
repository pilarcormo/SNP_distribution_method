require 'test/unit'
require_relative '../lib/reform_ratio'
require_relative '../lib/output'

class TestOutput < Test::Unit::TestCase	
	def test_create_perm_fasta
		perm = []
		@fasta_array = ReformRatio::fasta_array("test/test2.fasta")
		ids, lengths = ReformRatio.fasta_id_n_lengths(@fasta_array)
		perm = ["frag1", "frag3", "frag2"]
		fasta_perm = Output.create_perm_fasta(perm, @fasta_array, ids)
		assert_equal(fasta_perm, [">frag1 Length = 8", "CCAAATAC\n", ">frag3 Length = 7", "ACGACAC\n", ">frag2 Length = 8", "GCAATCGG\n"])
	end 
end 