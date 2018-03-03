import numpy as np

def compare(d):
	seqs=[d[k] for k in d]
	s1=seqs[0]
	s2=seqs[1]
	n=len(s1)
	m=len(s2)
	max_len, l, count,new_i,new_j=0,0,0,0,0
	bla,bla2,bla3,uniq_s2,match1,match2=[],[],[],[],[],[]
	score1=np.zeros(n)
	score2=np.zeros(m)
	for i,seq1 in enumerate(s1):
		for j,seq2 in enumerate(s2):
			#print new_j,j,l,max_len
			if seq1 == seq2  and (score1[i]==0 or score2[j]==0 ):
				
				count+=1
				l=len(match1)
				new_i=i+1
				new_j=j+1
				if score1[i]==0:
					match1.append(s1[i])
					score1[i]=1
					bla.append(seq1)
				if score1[j]==0:
					match2.append(s2[j])
					score2[j]=1
				if l > max_len:
					max_len=l
				print max_len
				print "match: score1[i] score2[j]"
				print score1[i],score2[j]
				print s1[i], s2[j]
				print " ---------- "
				break

			else:
				if new_i>=i:
					print "NOT A match: i,j,seq1,seq2"
					print score1[i],score2[j]
					print s1[i], s2[j]
					print " ---------- "
					if score2[j]==0:
						match2.append(s2[j])
						score2[j]=1
					else:
						match2.append("-")
				if new_j>=j:
					if score1[i]==0:
						match1.append(s1[i])
						score1[i]=1
					else:
						match1.append("-")
					bla2.append(seq1)
				l=0
				print max_len
				max_len=0
		uniq_s2.append(seq2)
	print bla
	print bla2
	print uniq_s2
	print match1
	print match2
	print max_len


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

def main():
	filename="bla"
	#with open(filename,"r") as f:
	compare(read_fasta(filename))

if __name__ == "__main__":
    main()