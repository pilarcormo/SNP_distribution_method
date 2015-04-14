#encoding: utf-8
require 'csv'


CSV.open("arabidopsis_datasets/#{dataset}/ratio_positions#{threshold}_#{adjust}.csv", "wb") do |csv|
  csv << ["Position", "Ratio"]
end

dic_poshm_short = {}
dic_poshm_short = dic_pos_hm

dic_poshm_short.each do |id, array|
  if dic_ratios.has_key?(id)
  else 
    dic_poshm_short.delete(id)
  end
end

dic_poshm_short.each do |id, array|
  array.each do |elem| 
    CSV.open("arabidopsis_datasets/#{dataset}/ratio_positions#{threshold}_#{adjust}", "ab") do |csv|
      csv << [elem, dic_ratios[id]] 
    end
  end 
end 