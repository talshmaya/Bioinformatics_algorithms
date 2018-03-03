proteins = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
    "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
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

def translate(DNA):
	seq=[]
	seqs=[]
	for index in range(0,len(str(DNA))):
		if DNA[index:index+3] == "AUG":	
			if len(seq)>0:
				#print ''.join(seq)
				seqs.append(''.join(seq))
			i=index		
			seq=[]
			for indexx in range(i,len(str(DNA)),3):	
				if len(DNA[indexx:indexx+3]) == 3 and proteins[DNA[indexx:indexx+3]] != "":
					seq.append(proteins[DNA[indexx:indexx+3]])
					
				else:
					break

	return seqs

def complementry(line):
	new_line = line.replace('A', 'u').replace('U', 'a').replace('G', 'c').replace('C', 'g')
	return new_line.upper()

def main():
	#with open ("rosalind_orf.txt","r") as f:
#		next(f)				
#		for line in f:
			#complementry(line)
	for k,line in read_fasta("rosalind_orf.txt").iteritems():
		#print line
		line = line.replace('T', 'U')
		line.strip()
		forward=translate(line.strip())
		reverse=translate(complementry(line[::-1].strip()))
		#print reverse
	print "\n".join(set(forward+reverse))
if __name__ == "__main__":
    main()
