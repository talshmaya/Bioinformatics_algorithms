class Node:
	""" Node class bla bla """
	def __init__(self,kmer):
		self.kmer=kmer
		self.edges=set()

	def get_kmer(self) :
		return self.kmer

	def print_me(self ):
		arr = []
		for i in self.edges:
			arr.append(i.get_kmer())
		print(self.kmer, arr)
	def __str__(self):
		return str(self.kmer)
		#self.kmer=kmer
		#return "kmer={0}/nedges={1}".format(self.kmer,self.edges)

def complementry(line):
	new_line = line.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
	return new_line.upper()

def split2kmers(seq,kmer):
		l=len(seq)
		b=[]
		#for i in range(0,l-kmer+1):
		#	b.append(seq[i:i+kmer])
		b.append(seq[:l-1])
		b.append(seq[1:])
		return b
def main():
	output_file= open('dbruj_out.txt', 'r+')
	with open ("rosalind_dbru.txt","r") as f:
	#with open ("bla.txt","r") as f:
		node_dict={}
		
		lines=set()
		for line in f:
			line=line.strip()
			km=len(line)-1
			c=complementry(line[::-1])
			lines.add(line)
			lines.add(c)
			
	for reads in lines:
		last_node=None
		for k in split2kmers(reads,km):
			if k not in node_dict:
				node_dict[k]=Node(k)
			if last_node:
				last_node.edges.add(node_dict[k])
			last_node=node_dict[k]

	for k,v in sorted(node_dict.iteritems()):
		for edge in sorted(v.edges):
			#if str(k)!=str(edge):
			output_file.write("({0}, {1})".format(k,edge)+'\n')
			print "({0}, {1})".format(k,edge)
	output_file.close()

if __name__ == "__main__":
    main()