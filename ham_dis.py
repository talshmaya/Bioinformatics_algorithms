def compare(s1,s2):
	count=0
	for s1,s2 in zip(s1,s2):
		if s1!=s2:
			count+=1
	return count

def main():
	filename="rosalind_hamm.txt"
	with open(filename,"r") as f:
		lines=f.readlines()
		seq1=lines[0]
		seq2=lines[1]
	print compare(seq1,seq2)

if __name__ == "__main__":
    main()
