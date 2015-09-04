
for i in {1..5} 
do
	ruby remove_cent.rb $i Aw_sup1-2/filter2_chromosome$i
	ruby remove_cent.rb $i BCF2/BCF2_chromosome$i
	ruby remove_cent.rb $i OCF2/OCF2_chromosome$i
	ruby remove_cent.rb $i B/B_chromosome$i
	ruby remove_cent.rb $i C/C_chromosome$i
done
