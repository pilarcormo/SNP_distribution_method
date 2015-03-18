

desc "SDM"
task :SDM  do
	sh "ruby SNP_distribution_method_variation.rb 'BCF2_2' 'blabla'"
end

file 'lolo.txt' => ["SDM"] do
	sh 'touch lolo.txt'
end 


desc "Convert sam to bam file"
task :samtools  do
	sh "arabidopsis_datasets/BCF2_4/blabla/mutation.txt"
end



# ['SDM', 'samtools'].each do |task|
#   Rake::Task[task].invoke
# end

