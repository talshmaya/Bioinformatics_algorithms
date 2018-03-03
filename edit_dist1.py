import numpy as np

def compare(d):
	seqs=[d[k] for k in d]
	s1="-"+str(seqs[0])
	s2="-"+str(seqs[1])
	n=len(s1)
	m=len(s2)
	score,max_score=0,0
	score_mx=np.zeros((m,n))
	traceback_mx=np.zeros((m,n))
	for i,seq1 in enumerate(s2):
			for j,seq2 in enumerate(s1):
					if i!=0 and j!=0:
						if s1[j]==s2[i]: #match
							score=1
						else: # mismatch 
							score=0
						
						score_mx[i][j]=max(score_mx[i-1][j-1]+score,score_mx[i-1][j],score_mx[i][j-1])
						
						if score_mx[i][j]>max_score:
							max_score=score_mx[i][j]
						if max(score_mx[i-1][j-1]+score,score_mx[i-1][j],score_mx[i][j-1],0)==score_mx[i-1][j]:
							traceback_mx[i][j]=2
						if max(score_mx[i-1][j-1]+score,score_mx[i-1][j],score_mx[i][j-1],0)==score_mx[i][j-1]:
							traceback_mx[i][j]=1
						
						if max(score_mx[i-1][j-1]+score,score_mx[i-1][j],score_mx[i][j-1],0)==score_mx[i-1][j-1]+score:
							traceback_mx[i][j]=3
	new_s1,new_s2=[],[]
	new_i,new_j=m-1,n-1
	for i in reversed(range(1,m)):
		
		for j in reversed(range(1,n)):
			if i==new_i and j==new_j:

				if traceback_mx[i][j]==2:
					new_i=i-1
					new_j=j
					new_s1.append("-")
					new_s2.append(s2[i])
				if traceback_mx[i][j]==1:
					new_i=i
					new_j=j-1 
					new_s1.append(s1[j])
					new_s2.append("-")
				if traceback_mx[i][j]==3:
					new_i=i-1
					new_j=j-1
					new_s1.append(s1[j])
					new_s2.append(s2[i])
	
	print int(max_score)
	print ''.join(new_s1[::-1])
	print ''.join(new_s2[::-1])
	print traceback_mx

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
	#filename="rosalind_edit.txt"
	filename="Q17.txt"
	compare(read_fasta(filename))

if __name__ == "__main__":
    main()