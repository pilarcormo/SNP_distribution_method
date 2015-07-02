#encoding: utf-8
class PDist

	require 'diff/lcs'

	# Returns array of values for each object's deviation in permutation, from original
	def self.distances(original, permutation)
		#################Change################
		indices = []
		indices_original = []
		original.each do |i|
		if permutation.include?(i) then
			indices_original << permutation.index(i) # indices of original values in permutation
		else 
			indices_original << original.length
		end
		indices << indices_original
		indices << Array(0..(original.length - 1))
		difference = indices.map {|x| x.reduce(:-)}
		return difference.map! {|i| i.abs } # the deviation's, as absolute values (direction does not matter)
		#####################################

		############Previous###################
		# indices = []
		# indices << original.map{|x| permutation.index(x)} # indices of original values in permutation
		# indices << Array(0..(original.length - 1)) # indices of original values in original
		# difference = indices.transpose.map {|x| x.reduce(:-)} # subtracting each object's index in original, from it's index in permutation, to find deviation
		# return difference.map! {|i| i.abs } # the deviation's, as absolute values (direction does not matter)
		#######################################
	end

	# Returns float of the deviation distance between original and permutation
	def self.deviation_raw(original, permutation)
		return distances(original, permutation).inject(:+) # sum of each object's deviation, is the deviation distance between original and permutation
	end

	# Returns float of the deviation distance between original and permutation, normalized between 0.0 and 1.0
	def self.deviation(original, permutation)
		s = deviation_raw(original, permutation)
		n = permutation.length
		if n % 2 == 0 
			distance = (2.0 / (n ** 2).to_f) * s
		else
			distance = (2.0 / ((n ** 2) - 1).to_f) * s
		end
		return distance # if/else block normalizes the deviation distance
	end
end
# 	# Returns float of the squared deviation distance between original and permutation
# 	def self.square_raw(original, permutation)
# 		sq_dists = []
# 		distances(original, permutation).each{|d| sq_dists << d**2} # each object's deviation is squared
# 		return sq_dists.inject(:+) # sum of each object's squared deviation, is the squared deviation distance between original and permutation
# 	end

# 	# Returns float of the squared deviation distance between original and permutation, normalized between 0.0 and 1.0
# 	def self.square(original, permutation)
# 		s = square_raw(original, permutation)
# 		n = permutation.length
# 		return (3.0 / (n**3 - n).to_f) * s # normalizes the squared deviation distance
# 	end

# 	# Returns float of the generalized hamming distance between original and permutation
# 	def self.hamming_raw(original, permutation)
# 		x = 0
# 		hds = [] # hamming distances
# 		permutation.each do |object| # hamming distances are 0 when object's have same index in original and permutation, and 1 when not
# 			if object != original[x]
# 				hds << 1
# 			end
# 			x+=1
# 		end
# 		return hds.inject(:+).to_f # total hamming distance is the sum of the hamming distances
# 	end

# 	# Returns float of the generalized hamming distance between original and permutation, normalized between 0.0 and 1.0 
# 	def self.hamming(original, permutation)
# 		total_hd = hamming_raw(original, permutation)
# 		return total_hd / permutation.length.to_f # normalizes by dividing the total hamming distance by the maximum possible total, which == number of objects
# 	end

# 	# Returns float of the compliment R distance between original and permutation
# 	def self.rdist_raw(original, permutation)
# 		x = 0
# 		r = [] # compliment of R distance == number of times two objects consecutive in original, are consectutive in permutation
# 		(permutation.length - 1).times do 
# 			y = permutation.index(original[x])
# 			if original[x+1] != permutation[y+1]
# 				r << 1
# 			end
# 			x+=1
# 		end
# 		return r.inject(:+).to_f
# 	end

# 	# Returns float of the compliment R distance between original and permutation, normalized between 0.0 and 1.0
# 	def self.rdist(original, permutation)
# 		dist = rdist_raw(original, permutation)
# 		return dist / (permutation.length - 1).to_f # normalized: dividing by maximum number of consecutive occurences
# 	end

# 	# Returns float of the longest common sub-sequence between original and permutation
# 	def self.lcs_raw(original, permutation)
# 		return Diff::LCS.LCS(original, permutation) # diff-lcs gem used to calculate longest common sub-sequence
# 	end

# 	# Returns float of the longest common sub-sequence between original and permutation, normalized between 0.0 and 1.0
# 	def self.lcs(original, permutation)
# 		lcs = lcs_raw(original, permutation)	
# 		return (permutation.length - lcs.length).to_f / (permutation.length - 1).to_f # normalized: dividing by longest possible common sub-sequence
# 	end

# 	# Returns float of the kendall's tau distance between original and permutation
# 	def self.kendalls_tau_raw(original, permutation)
# 		kt = []
# 		original.each do |x| # for each of the objects in original...
# 			permutation.each do |y| # ... iterate over the objects in permutation
# 				# every time object y comes before object x in permutation, where it comes after in original... 
# 				if original.index(x) < original.index(y) && permutation.index(x) > permutation.index(y)
# 					kt << 1 # ... 1 is added
# 				end
# 			end
# 		end
# 		return kt.inject(:+) # the number of pairwise adjacent permutations required to transform original into permutation
# 	end

# 	# Returns float of the kendall's tau distance between original and permutation, normalized between 0.0 and 1.0
# 	def self.kendalls_tau(original, permutation) # TODO edit the if statement to work with any type of object
# 		n = permutation.length
# 		s = kendalls_tau_raw(original, permutation)
# 		return 2 * (s.to_f / (n**2 - n).to_f) # normalized kendall's tau distance
# 	end
end