while read line
do
	#bla=$line-2
	head -$line pza_random_mutations2 |tail -3
done < $1