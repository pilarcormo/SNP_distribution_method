#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require 'pp'
require "csv"


dataset = ARGV[0]

fasta = "arabidopsis_datasets/#{dataset}/frags_ordered1.fasta"

frags = ReformRatio.fasta_array(fasta)

ids_ok, lengths_ok, id_len_ok = ReformRatio.fasta_id_n_lengths(frags)

pp ids_ok

pos = []
csv = "arabidopsis_datasets/#{dataset}/ratio_positions0_1.csv"


dic = {}
density_pos = []

CSV.foreach(csv, :headers => true, :header_converters => :symbol, :converters => :all) do |row|
  dic.store(row.fields[0],row.fields[1])
end

dic.each do |pos, ratio|
	(ratio*10).to_i.times do
		density_pos << pos 
	end 
end 

WriteIt::write_txt("arabidopsis_datasets/#{dataset}/ratios", density_pos)