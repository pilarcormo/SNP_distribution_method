directory "Pilar"

desc "Print hello"
file "hello.txt" => "Pilar" do 
	sh "echo 'hello' >> 'Pilar/hello.txt'"
end

task :test do 
	sh "ls"
end

# desc "builds output sam"
# file "new_BCF2.sam" => "BCF2" do 
# 	sh "bwa.sh"
# end 

# desc "runs blabla"
# file 'hello.txt' => 'output.sam' do 
# 	sh 'echo "running samtools"'
# end 

# desc "makes vcf from sam"
# file "vcf.vcf "

# task :default => "vcf.vcf" 

desc "sam"
file "output.sam" => "hello.txt" do 
	sh 'echo "hello"'
end

task :sam => "output.sam" do 
	puts 'blabla'
end 



task :caca => [:sam] do 
	puts 'lele'
end 


