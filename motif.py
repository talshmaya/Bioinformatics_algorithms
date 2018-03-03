def compare(s,m):
	count=0
	hit=[]
	m_len=len(m) # length of the motif
	s_len=len(s)
	for i in range(0,s_len-m_len+1):
		if str(s[i:i+m_len-1]).strip() == str(m).strip():
			hit.append(str(i+1)) 
	return ' '.join(hit)

def main():
	filename="rosalind_subs.txt"
	with open(filename,"r") as f:
		lines=f.readlines()
		seq=lines[0]
		motif=lines[1]
	print compare(seq,motif)

if __name__ == "__main__":
    main()
