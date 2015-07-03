require_relative "lib/write_it"
require 'pp'
file = ARGV[0]


list = WriteIt.file_to_array("arabidopsis_datasets/#{file}/mutation.txt")

a = list.last
b = a.split(" ")
puts b[1]

File.open("arabidopsis_datasets/Genomes_SDM/30Mb/mutation.txt", "w+") do |f|
	f << b[1]
end

	