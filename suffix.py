def bwt(s):
	count=0
	seq_index=[]
	for i,char in enumerate(range(0,len(s)-1)):
		a=str(s[count:]).strip()
		b=str(s[:count]).strip()
		seq_index.append((a+b,i))
		count+=1
	indexes=[ str(i[1]) for i in sorted(seq_index)]
	print ', '.join(indexes)
	#for i in sorted(seq_index):
	#	print i[1]

def main():
	filename="bla.txt"
	with open(filename,"r") as f:
		lines=f.readlines()
		seq=lines[0]
	bwt(seq)

if __name__ == "__main__":
    main()

