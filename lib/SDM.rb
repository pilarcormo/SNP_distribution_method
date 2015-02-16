#encoding: utf-8

class SDM
	def self.sorting(dic_hm_inv)
		list1_hm, list2_hm, left_hm, right_hm = [], [], [], []
		keys_hm = dic_hm_inv.keys.to_a
		Array(1..keys_hm.length/2).each do |i|
      min1 = keys_hm.min
      list1_hm << dic_hm_inv.values_at(min1)
      list1_hm.flatten!
      keys_hm.delete(min1)
      if list1_hm.length.to_i % 2 == 0
        lu = list1_hm.each_slice(list1_hm.length/2.to_i).to_a
        right_hm << lu[0]
        left_hm << lu[1]
      else
        if list1_hm.length.to_i > 2
          object = list1_hm.shift
          lu2 = list1_hm.each_slice(list1_hm.length/2.to_i).to_a
          right_hm << lu2[0]
          left_hm << lu2[1]
          right_hm << object
        else
          right_hm << list1_hm
        end
      end
      min2 = keys_hm.min
      keys_hm.delete(min2)
      list2_hm << dic_hm_inv.values_at(min2)
      list2_hm.flatten!
      if list2_hm.length.to_i % 2 == 0
        lu = list2_hm.each_slice(list2_hm.length/2.to_i).to_a
        right_hm << lu[0]
        left_hm << lu[1]
      else
        if list2_hm.length.to_i > 2
          object = list2_hm.shift
          lu2 = list2_hm.each_slice(list2_hm.length/2.to_i).to_a
          right_hm << lu2[0]
          left_hm << lu2[1]
          left_hm << object
        else
          left_hm << list2_hm
        end
      end
      list1_hm = []
      list2_hm = []
    end
		right_hm = right_hm.flatten
		left_hm = left_hm.flatten.compact
		left_hm = left_hm.reverse #we need to reverse the left array to build the distribution properly
		perm_hm = right_hm << left_hm #combine together both sides of the distribution
		perm_hm.flatten!
		return perm_hm
	end
end