#!/usr/bin/perl
# generate an NGM file for mapping from VCF or BCF data as generated with 'samtools mpileup'
#
# Output is gzipped file with format:
# chr	pos	ref	alt	AF	DP	MQ	FQ	QUAL
#
# Discordant chastity is replaced with allele frequency as calculated from DP4 high quality counts
#
# Author: Ryan S. Austin
# Agriculture & AgriFood Canada
# Copyright January 2013, All Rights Reserved
#
$requires = 'bcftools from Samtools >v0.1.16, gzip';
$version = '0.1.1';  # beta

unless (@ARGV)
	{
	print STDERR "BCF2NGM v$version\nusage: bcf2ngm <file.bcf>\nBCF generated with: 'samtools mpileup -E -ugf <reference.fasta> <alignment.bam> | bcftools view -bvgN - > output.bcf'\nRequires: $requires\n";
	exit 1;
	}
$vcf = shift(@ARGV);
( -e "/tmp") ? $p = "/tmp" : $p = ".";
$tmp = "$p/tmp.vcf2ngm.$$";
$tmpZ = "$p/tmpZ.vcf2ngm.$$";
open(TMP,">$tmp");
if (`file $vcf | grep compressed`)
	{
	die "error: unable to find bcftools in your path, is samtools installed?\n" if system("which bcftools >/dev/null");
	system("bcftools view $vcf > $tmpZ");
	open(VCF, $tmpZ) or die "error: failed to open file $tmpZ : $!\n";
	}
else
	{
	open(VCF,$vcf) or die "error: failed to open file $vcf : $!\n";
	}
while(<VCF>)
	{
	next if $_ =~ /^#/;
	next if $_ =~ /INDEL/;   # this should be incorporated eventually 
	@f = split(/\s+/,$_);
	@info = split(/;/,$f[7]);
	$DP = $MQ = $FQ = $AF =  undef;
	foreach $x (@info)          # DP, MQ and FQ are as defined in samtools VCF header
		{
		if ($x =~ /DP=(\d+)/) # depth of high quality reads 
			{
			$DP = $1;
			}
		elsif ($x =~ /MQ=(-?\d+)/)  # mapping quality 
			{
			$MQ = $1;
			}
		elsif ($x =~ /FQ=(-?\d+)/)  # probability of all samples being the same (not used) 
			{
			$FQ = $1;
			}
		elsif ($x =~ /DP4=(\d+),(\d+),(\d+),(\d+)/)  # high quality refForward, refReverse, altForward, altReverse reads  
			{
			$AF = ($3 + $4) / ($1 + $2 + $3 + $4);  # reference allele frequencey - replaces chastity 
			$AF = sprintf("%1.2f",$AF);
			}
		else
			{
			next;
			}
		}
	$QUAL = $f[5];
	die "error: failed to parse DP=$DP/MQ=$MQ/FQ=$FQ/AF=$AF from:\n $_" unless ($DP and defined $FQ and defined $MQ and defined $AF);
	$chromosome = $f[0];
	$chromosome =~ s/chromosome\s*//i;
	$chromosome =~ s/chr\s*//i;
	$chromosome =~ s/scaffold\s*//i;
	$position = $f[1];
	$referenceBase = $f[3];
	$altBase = $f[4];
	unless ($referenceBase eq 'N')
		{
		print TMP "$chromosome\t$position\t$referenceBase\t$altBase\t$AF\t$DP\t$MQ\t$MQ\t$QUAL\n";
		}
	}
system("gzip -c $tmp > $vcf.ngm");
unlink($tmp);
unlink($tmpZ);
print "Finished.  Output written to $vcf.ngm\n";
exit 0;
