import numpy as np
from collections import OrderedDict

def read_fasta(filename):
	my_dict=OrderedDict()
	header=""
	with open(filename,"r") as f:
		lines=f.readlines()
	for i,line in enumerate(lines):
		if ">" in line:
			header=line.strip()
		else:
			my_dict[header]=str(my_dict.get(header,""))+str(line.strip())
			
	return my_dict	

def compare(seqs,matrix):

	for i,seq1 in enumerate(seqs):
		for j,seq2 in enumerate(seqs):
			score=distance(seq1,seq2)/float(len(seq1))
			matrix[i][j]='{0:.5f}'.format(float(score))

	return matrix

def distance(s1,s2):
	count=0
	for ss1,ss2 in zip(s1,s2):
		if ss1!=ss2:
			count+=1
	return count

def main():
	output_file='rosalind_pdst.txt'
	d=read_fasta(output_file)

	num_of_seqs=len(d)
	prev_seq=[]
	score_mx=np.zeros((num_of_seqs,num_of_seqs),dtype='object')

	for i in compare(d.values(),score_mx):
		print '\t'.join(i) 

if __name__ == "__main__":
    main()
