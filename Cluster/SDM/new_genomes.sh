source ruby-2.0.0
source R-3.1.0



for i in {1..5} 
do
	bsub -q TSL-Test128 -o 1.txt "ruby model_genome_new.rb 1Mb_$i 1000000 1000"
	bsub -q TSL-Test128 -o 1A.txt "ruby model_genome_new.rb 1Mb_A_$i 1000000 500"
	bsub -q TSL-Test128 -o 3.txt "ruby model_genome_new.rb 3Mb_$i 3000000 3000"
	bsub -q TSL-Test128 -o 3A.txt "ruby model_genome_new.rb 3Mb_A_$i 3000000 1500"
	bsub -q TSL-Test128 -o 3.txt "ruby model_genome_new.rb 5Mb_$i 5000000 5000"
	bsub -q TSL-Test128 -o 3A.txt "ruby model_genome_new.rb 5Mb_A_$i 5000000 2500"
	bsub -q TSL-Test128 -o 7.txt "ruby model_genome_new.rb 7Mb_$i 7000000 7000"
	bsub -q TSL-Test128 -o 7A.txt "ruby model_genome_new.rb 7Mb_A_$i 7000000 3500"
	bsub -q TSL-Test128 -o 9.txt "ruby model_genome_new.rb 9Mb_$i 9000000 9000"
	bsub -q TSL-Test128 -o 9A.txt "ruby model_genome_new.rb 9Mb_A_$i 9000000 4500"
	bsub -q TSL-Test128 -o 11.txt "ruby model_genome_new.rb 11Mb_$i 11000000 11000"
	bsub -q TSL-Test128 -o 11A.txt "ruby model_genome_new.rb 11Mb_A_$i 11000000 5500"
	bsub -q TSL-Test128 -o 13.txt "ruby model_genome_new.rb 13Mb_$i 13000000 13000"
	bsub -q TSL-Test128 -o 13A.txt "ruby model_genome_new.rb 13Mb_A_$i 13000000 6500"
	bsub -q TSL-Test128 -o 15.txt "ruby model_genome_new.rb 15Mb_$i 15000000 15000"
	bsub -q TSL-Test128 -o 15A.txt "ruby model_genome_new.rb 15Mb_A_$i 15000000 7500"
done 
