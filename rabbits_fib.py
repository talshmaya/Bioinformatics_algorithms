def rabbits(m,n):
	rabbits=[1]
	for i in range(0,m-1):
		rabbits.append(0)
	for i in range(0,n-1):
		rabbits=[sum(rabbits[1:])]+rabbits[:-1]
	return sum(rabbits)

print rabbits(20,83)