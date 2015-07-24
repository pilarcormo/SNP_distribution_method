require_relative "lib/write_it"
require 'pp'
file = ARGV[0]

mutation = WriteIt.file_to_array("arabidopsis_datasets/#{file}/mutation.txt")

deviation = mutation.last
percentage = deviation.split(" ")
puts percentage[1]

File.open("arabidopsis_datasets/Genomes_SDM/30Mb/mutation.txt", "w+") do |f|
	f << percentage[1]
end

	