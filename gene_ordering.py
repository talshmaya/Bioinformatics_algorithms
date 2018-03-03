from itertools import *
from collections import defaultdict

# def perm(n):
# 	positives=list(permutations(range(1,n+1)))
# 	negatives=list(permutations(range(-1,-n-1,-1)))
# 	counter,passed=[],[]
# 	for r in product(positives, negatives): 
# 		for i in range(0,len(r)):
# 			first_pos=r[0][:i+1]
# 			second_neg=r[1][i+1:]
# 			first_neg=r[1][:i+1]
# 			second_pos=r[0][i+1:]
# 			l=first_pos+second_neg
# 			counter.append(' '.join(str(char) for char in first_pos+second_neg))
# 			counter.append(' '.join(str(char) for char in first_neg+second_pos))
# 			sett=set(counter)
# 	for k in sett:
# 		digs=set()
# 		for m in k:
# 			if m.isdigit():
# 				digs.add(abs(int(m)))
# 		if len(digs)==n:
# 			passed.append(k)

# 	print len(passed)
# 	for i in passed:
# 		print i

def perm(n):
	my_list = range(1,n+1) + range(-n,0)
	arr = []
	for v in permutations(my_list, n):
		d = defaultdict(int)
		for j in v:
			if j in d:
				pass
			else:
				d[abs(j)] += 1

		if len(d.keys()) == n:
			arr.append(v)
	print(len(set(arr)))
	for i in arr:
		print(str(i).replace('(','').replace(')', '').replace(',', ''))
		# print(d)


def main():
	perm(5)

if __name__ == "__main__":
    main()
