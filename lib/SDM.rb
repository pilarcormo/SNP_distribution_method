#encoding: utf-8


class SDM
	def self.sorting(dic_hm_inv)
    def self.divide_array(dic_hm_inv, right, left, keys_hm, dest)
      contigs_at_min = []
      minimum = keys_hm.min
      contigs_at_min << dic_hm_inv.values_at(minimum)
      contigs_at_min.flatten!
      keys_hm.delete(minimum)
      if contigs_at_min.length.to_i % 2 == 0
        half = contigs_at_min.each_slice(contigs_at_min.length/2.to_i).to_a
        right << half[0]
        left << half[1]
      else
        if contigs_at_min.length.to_i > 2
          object = contigs_at_min.shift
          other_half = contigs_at_min.each_slice(contigs_at_min.length/2.to_i).to_a
          right << other_half[0]
          left << other_half[1]
          right << object
        else 
          if dest == 0
            right << contigs_at_min
          elsif dest == 1
            left << contigs_at_min
          end
        end
      end
      return right, left, keys_hm
    end
    left, right = [], []  
		keys= dic_hm_inv.keys.to_a
		Array(1..keys.length/2).each do |i|
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 0)
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 1)
    end 
    perm = right.flatten << left.compact.flatten.reverse #combine together both sides of the distribution
    perm.flatten!
    mut = []
    mut << right.flatten[-2, 2]
    mut << left.flatten[-2, 2].reverse
    mut.flatten!
    return perm, mut
	end
end