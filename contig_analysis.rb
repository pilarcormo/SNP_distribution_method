#encoding: utf-8
require 'pp'
require_relative 'lib/write_it'
require 'rinruby'
require 'csv'
#N50 = "Given a set of contigs, each with its own length, the N50 length is defined as the length for which the collection of all contigs of that length or longer contains at least half of the sum of the lengths of all contigs, and for which the collection of all contigs of that length or shorter also contains at least half of the sum of the lengths of all contigs. "

csv_array = WriteIt.file_to_array("Contigs/contigs.txt")

species, lengths, contigs, n_length, coverage, technology = [], [], [], [], [], []
csv_array.each do |line|
	a = line.split("|")
	species << a[0] 
	lengths << a[1]
	contigs << a[2]
	n_length << a[3]
	coverage << a[4]
	technology << a[5]
end 

species.shift 
lengths.shift
contigs.shift
n_length.shift
coverage.shift
technology.shift

average_contig, average_contig_n50 = [], []

Array(0..lengths.length-1).each do |x|
	average_contig << lengths[x].to_i/contigs[x].to_i
	average_contig_n50 << lengths[x].to_i/n_length[x].to_i
end 

CSV.open("contigs.csv", "wb") do |csv|
	csv << ["Species", "Sequence_lengths", "Number_of_contigs", "N50_contig", "Average_contig_size", "Average_N50", "Coverage", "Technology_used"]
end 

Array(0..lengths.length-1).each do |x|
	CSV.open("contigs.csv", "ab") do |csv|
		csv << [species[x], lengths[x], contigs[x], n_length[x], average_contig[x], average_contig_n50[x], coverage[x], technology[x]]
	end 
end 


