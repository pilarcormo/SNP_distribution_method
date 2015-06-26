#encoding: utf-8


class SDM
  ##Input for sorting: inverted hash containing the normalised homozygous scores as key and the fragments' ids as value and the type of cross (back or out). 
  ##Output from sorting: perm is the permutation array containing the candidate order of contigs and mut is the array of candidate contigs taken from the central part of perm
	def self.sorting(dic_hm_inv, cross)
    def self.divide_array(dic_hm_inv, right, left, keys_hm, dest)
      contigs_at_min = []
      minimum = keys_hm.min #define minimum score 
      contigs_at_min << dic_hm_inv.values_at(minimum) #look for the contigs with the minimum score 
      contigs_at_min.flatten!
      keys_hm.delete(minimum) #delete entries for contigs with the minimum homozygous score from the original hash
      #shuffle the contigs with the minimum homozygous scores on the two halfs of the expected normal distribution. 
      #(1) If the number of contigs is even, half of the array goes to the left and half goes to the right part of the distribution
      #(2) If the number of contigs is odd, we define two situations:
            ##if there's only 1 contig, we randomly assign the right or left destination for it
            ##if there are more than 1 contigs, we take the first contig in the array and randomly assign the right or left destination for it. The remaining array have a even number of elements, so we proceed as described in (1)
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
		Array(1..keys.length/2).each do |i| #repeat the sorting process until the original hash is sorted. 
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 0)
      right, left, keys = SDM.divide_array(dic_hm_inv, right, left, keys, 1)
    end 
    perm = right.flatten << left.compact.flatten.reverse #combine together both sides of the distribution
    perm.flatten!
    mut = []
    pl = perm.length
    pp right
    pp left
    if cross == "back" #In case of a backcross, 10 contigs in the middle part of the permutation taken
      mut << right.flatten[-5, 5]
      mut << left.flatten[-5, 5].reverse
      mut.flatten!
    elsif cross == "out" #In case of a backcross, 20 contigs in the middle part of the permutation taken
      if pl > 20
        mut << right.flatten[-10, 10]
        mut << left.flatten[-10, 10].reverse
        mut.flatten!
      else #If a strong filtering step reduces the total number of contigs to a number lower than 20, perm.length/2 contigs on the right and perm.length/2 on the left side of the middle point are taken. 
        mut << right.flatten[-pl/2, pl/2]
        mut << left.flatten[-pl/2, pl/2].reverse
        mut.flatten!
      end 

    end 
    pp mut 
    return perm, mut
	end

  ##Input for calling: 
  ##Output from calling:
  def self.calling_SDM (dic_hm_inv, dic_ratios_inv, cross, dic_pos_hm)
    mut, number_of_snps = [], []
    perm_hm, mut  = SDM.sorting(dic_hm_inv, cross)
    perm_ratio, mut_ratio = SDM.sorting(dic_ratios_inv, cross)
    mut << mut_ratio
    mut.flatten!
    mut.uniq 
    or_pos = {}
    mut.each { |frag|
        if dic_pos_hm.has_key?(frag)
          or_pos.store(frag, dic_pos_hm[frag])
        end
      }
    or_pos.each do |frag, array|
        number_of_snps << array.length
    end
    #or_pos.delete_if { |id, array|  array.length < (number_of_snps.max - 1) }  
    hyp = or_pos.values.sort.flatten!

    return perm_hm, perm_ratio, mut, hyp 
  end 
end