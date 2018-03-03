import re
from collections import OrderedDict

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

def translate(DNA):
	seq=[]

	for i in range(0,len(str(DNA)),3):	
		if len(DNA[i:i+3]) == 3 and proteins[DNA[i:i+3]] != "":
			seq.append(proteins[DNA[i:i+3]])
			#print proteins[DNA[i:i+3]],DNA[i:i+3]

	return ''.join(seq)

def complementry(line):
	new_line = line.replace('A', 'u').replace('U', 'a').replace('G', 'c').replace('C', 'g')
	return new_line.upper()

def splice(introns,s):
	windows=[]
	for i in introns:
		#print len(i),len(s)
		while i in s:
			s=re.sub(i,"",s)
	return s

def main():
	introns=[]
	d=read_fasta("rosalind_splc.txt")
	seq=d.values()[0]
	key_seq=d.keys()[0]
	del d[key_seq]
	for k,line in d.iteritems():
		introns.append(line)
	spliced=splice(introns,seq)
	spliced = spliced.replace('T', 'U')
	spliced.strip()
	print translate(spliced.strip())
	
if __name__ == "__main__":
    main()

