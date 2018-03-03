def main():
	input_file='rosalind_ba10b.txt'
	with open(input_file,"r") as f:
		lines=f.readlines()
	
	seq=lines[0]
	ems=lines[4]

	A2x=float(lines[9].split()[1])
	A2y=float(lines[9].split()[2])
	A2z=float(lines[9].split()[3])

	B2x=float(lines[10].split()[1])
	B2y=float(lines[10].split()[2])
	B2z=float(lines[10].split()[3])

	#d = {
	#	'A' : {'x': A2x, 'y': A2y,'z': A2z},
	#	'B' : {'x': B2x, 'y': B2y,'z': B2z}
	#}

	prob=[]

	for i,j in zip(seq,ems):
		trans=j+i
		if trans=="Ax":
			prob.append(A2x)
		if trans=="Ay":
			prob.append(A2y)
		if trans=="Az":
			prob.append(A2z)
		if trans=="Bx":
			prob.append(B2x)
		if trans=="By":
			prob.append(B2y)
		if trans=="Bz":
			prob.append(B2z)
	
	print reduce(lambda x, y: x*y, prob)
if __name__ == "__main__":
    main()

