from itertools import permutations

def perm(n):
	#l=permutations(range(3))
	#for i in l:
	#	yield str(i)
	l=list(permutations(range(1,n+1)))
	print len(l)
	for i in l:
		print ' '.join(str(char) for char in i)

def main():
	filename="rosalind_edta.txt"
	perm(6)

if __name__ == "__main__":
    main()
