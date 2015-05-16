#encoding: utf-8
require_relative '../lib/vcf'
require_relative '../lib/write_it'
require 'test/unit'

class TestVCF < Test::Unit::TestCase
	def setup
		@vcf_ngs = "test/ngs.vcf"
		@chromosome = 1
		@vcfs_info = {"ADP"=>"17", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}, {"ADP"=>"25", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}
		@vcfs_pos = [5, 123]
		@snps = {5 => "HET", 123 => "HET"}
		@vcf = ["1\t5\t.\tC\tA\t.\tPASS\tADP=17;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:24:17:17:9:7:41.18%:3.3988E-3:65:52:9:0:1:6\n", 
		"1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"]
	end
	def test_open_vcf
		vcf, chrom, pos, info = Vcf.open_vcf(@vcf_ngs, @chromosome)
		assert_equal( ["1\t5\t.\tC\tA\t.\tPASS\tADP=17;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:24:17:17:9:7:41.18%:3.3988E-3:65:52:9:0:1:6\n", 
		"1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"], vcf)
		assert_equal(["1", "1"], chrom)
		assert_equal([5, 123], pos)
		assert_equal( [{"ADP"=>"17", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}, {"ADP"=>"25", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}], info)
	end 
	def test_type_per_pos
		snps, hm, ht = Vcf.type_per_pos(@vcfs_info, @vcfs_pos)
		assert_equal({5 => "HET", 123 => "HET"}, snps)
		assert_equal([], hm)
		assert_equal([5, 123], ht)
	end 
	def test_filtering
		snps_p = {5 => "HET", 365 => "HOM"}
		short_vcf = Vcf.filtering(@vcfs_pos, snps_p, @snps, @vcf)
		assert_equal(["1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"], short_vcf)
	end 
end 