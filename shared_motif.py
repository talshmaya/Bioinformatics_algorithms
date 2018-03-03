def read_fasta(filename):
	my_dict={}
	header=""
	with open(filename,"r") as f:
		lines=f.readlines()
	for i,line in enumerate(lines):
		if ">" in line:
			header=line.strip()
		else:
			my_dict[header]=str(my_dict.get(header,""))+str(line.strip())
			
	return my_dict

def split2kmers(seqs,kmer):
		
		b,bla=[],[]
		
		for seq in seqs:
			l=len(seq)
			b=[]
			#
			for i in range(0,l-kmer+1):
				b.append(seq[i:i+kmer])
			bla.append(b)
		return bla

def compare_seqs(seqs,num_of_seqs):
	common,prev_common="",""
	b,prev,c=[],[],[]
	for seq in seqs:
		common=set(seq).intersection(prev)
		if len(common)>0:
			c.append(common)
		prev=seq
		if len(c)==num_of_seqs-1:
			l=min(c)
			if l!="":
				b=l.pop()
	return b		

def main():
	output_file='bla.txt'
	d=read_fasta(output_file)
	num_of_seqs=len(d)
	prev_seq=[]
	max_len_seq=0
	for i in d.values():
		curr_len=len(i)
		if curr_len>max_len_seq:
			max_len_seq=curr_len
	
	for kmer in reversed(range(1,max_len_seq)):
		if len(compare_seqs(split2kmers(d.values(),kmer),num_of_seqs))>0:
			print compare_seqs(split2kmers(d.values(),kmer),num_of_seqs)
			break

if __name__ == "__main__":
    main()
