#encoding: utf-8

##Input1: permutation array after SDM
##Input2: fasta file converted to array
##Input3: list of fragment ids from the shuffled fasta file
##Output: permutation of fragments after SDM with the data (lengths, etc) obtained from the original fasta file
class Output 
	def self.create_perm_fasta(perm, frags, ids)
		defs, data, defs_p, data_p, fasta_perm = [], [], [], [], []
		frags.each do |i|
			defs << i.definition
			data << i.data
		end
		perm.each { |frag|
	      index_frag = ids.index(frag).to_i
	      defs_p << defs[index_frag]
	      data_p << data[index_frag]
    	}
    ###Create fasta array with the information above
		x = 0
		Array(0..perm.length-1).each do |i|
			fasta_perm << defs_p[x]
			fasta_perm << data_p[x]
			x += 1
		end
		fasta_perm.each do |line|
			if line.start_with?("\n")
				line[0] = ''
			else
				line.insert(0, ">")
			end
		end
		return fasta_perm
	end
end 
