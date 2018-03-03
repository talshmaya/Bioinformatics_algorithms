def bwt(s):
	count=0
	sorted_s=sorted(s.strip())
	while len(sorted_s[0])< len(s)-1:
		b=[i for i in sorted_s]
		sorted_s=[]
		a=[i for i in s.strip()]
		for aa,bb in zip(a,b):
			sorted_s.append(aa+bb)
		sorted_s=sorted(sorted_s)
	for i in sorted_s:
		if i[-1] == '$':
			return i

def main():
	filename="bla.txt"
	with open(filename,"r") as f:
		lines=f.readlines()
		seq=lines[0]
	print bwt(seq)

if __name__ == "__main__":
    main()

