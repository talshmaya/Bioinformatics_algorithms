def main():
	input_file='rosalind_ba10a.txt'
	with open(input_file,"r") as f:
		lines=f.readlines()
	
	seq=lines[0]
	states=lines[2]

	A2A=float(lines[5].split()[1])
	B2A=float(lines[6].split()[1])
	A2B=float(lines[5].split()[2])
	B2B=float(lines[6].split()[2])
	#print A2A, A2B, B2A, B2B
	
	prob=[]
	prob.append(0.5)
	
	for i in range(0,len(seq)-1):
		trans=seq[i:i+2]

		if trans=="AA":
			prob.append(A2A)
		if trans=="AB":
			prob.append(A2B)
		if trans=="BA":
			prob.append(B2A)
		if trans=="BB":
			prob.append(B2B)
	
	print reduce(lambda x, y: x*y, prob)
if __name__ == "__main__":
    main()

