#encoding: utf-8
require_relative 'stuff'
require 'csv'

class Ratio_filtering 

	def self.important_ratios(snps_hm, snps_ht, ids, threshold, adjust) 
		x = 0
		dic_ratios, ratios, ids_s = {}, [], []
		dic_ratios_by_length = {}
		snps_hm.length.times do
			ratio = (snps_hm[x]+adjust.to_f)/(snps_ht[x]+adjust.to_f)
			dic_ratios.store(ids[x], ratio.to_f) 
			x += 1
		end
		if threshold > 0
			thres = 100/threshold
			pp thres 
			filter = (dic_ratios.values.max.to_f)/thres
			puts "All the contigs with a ratio value below #{filter} will be discarded"
			dic_ratios.delete_if { |id, ratio|  ratio <= filter.to_f}
		end 
		ratios << dic_ratios.values
		ratios.flatten!
		ids_s = dic_ratios.keys
		dic_ratios_inv = Stuff.safe_invert(dic_ratios)
		return dic_ratios, ratios, ids_s, dic_ratios_inv
	end

	def self.important_ids(ids_short, ids)
		shuf_short_ids = []
		ids_short.flatten!
		ids.each do |frag|
			if ids_short.include?(frag)
		    	shuf_short_ids << frag
		  	end 
		end 
		shuf_short_ids.flatten!
		return shuf_short_ids
	end 
	def self.important_pos(ids_short, pos)
		sh = []
		pos.each do |frag, positions|
			if ids_short.include?(frag)
		  	else 
		    	pos.delete(frag)
		  	end 
		end 
		sh = pos.values
		sh.flatten!
		return sh
	end 

	#Input1 location for the csv file
	#Input2 hash with the id and positions for the hm SNPs
	#Input3 hash with the id and ratio for each fragment
	def self.csv_pos_ratio(csv, pos, ratios)
		pos_ratio = {}
		CSV.open(csv, "wb") do |csv|
		  csv << ["Position", "Ratio"]
		end
		short = pos
		short.each do |id, array|
		  if ratios.has_key?(id)
		  else 
		    short.delete(id)
		  end
		end
		short.each do |id, array|
			array.each do |elem| 
		    	CSV.open(csv, "ab") do |csv|
		      		csv << [elem, ratios[id]] 
		      		pos_ratio.store(elem, ratios[id])
		      	end 
		    end 
		end
		return pos_ratio 
	end 
