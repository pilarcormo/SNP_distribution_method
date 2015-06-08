#encoding: utf-8


class SDM
	def self.sorting(dic_hm_inv, cross)
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
    pl = perm.length
    if cross == "back"
      mut << right.flatten[-1, 1]
      mut << left.flatten[-1, 1].reverse
      mut.flatten!
    elsif cross == "out" 
      if pl > 10
        mut << right.flatten[-10, 10]
        mut << left.flatten[-10, 10].reverse
        mut.flatten!
      else 
        mut << right.flatten[-pl/2, pl/2]
        mut << left.flatten[-pl/2, pl/2].reverse
        mut.flatten!
      end 

    end 
    return perm, mut
	end
  
  def self.calling_SDM (dic_shuf_hm_norm, dic_ratios_inv_shuf, cross, dic_pos_hm)
    mut, number_of_snps = [], []
    perm_hm, mut  = SDM.sorting(dic_shuf_hm_norm, cross)
    perm_ratio, mut_ratio = SDM.sorting(dic_ratios_inv_shuf, cross)
    mut << mut_ratio
    mut.flatten!
    mut = mut.uniq!
    or_pos = {}
    mut.each { |frag|
        if dic_pos_hm.has_key?(frag)
          or_pos.store(frag, dic_pos_hm[frag])
        end
      }
    or_pos.each do |frag, array|
        number_of_snps << array.length
    end
    or_pos.delete_if { |id, array|  array.length < (number_of_snps.max - 1) }  
    hyp = or_pos.values.sort.flatten!
    return perm_hm, perm_ratio, mut, hyp 
  end 
end