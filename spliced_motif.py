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

def get_motif_index(seq,motif):
	
	new_i=0
	bla,mot=[],[]	
	for j,motif1 in enumerate(motif): # GTA
		for i,seq1 in enumerate(seq): # ACGTACGTGACG
			if i>new_i and seq1==motif1:
				new_i=i
				mot.append(motif1)
				bla.append(str(i+1))
				#print mot
				if ''.join(mot)==motif:
					print ' '.join(bla)
				break

def main():
	input_file='rosalind_sseq.txt'
	d=read_fasta(input_file)
	seq=d.values()[0]
	motif=d.values()[1]
	get_motif_index(seq,motif)

if __name__ == "__main__":
    main()

