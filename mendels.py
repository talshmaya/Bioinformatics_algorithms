def get_dominant(k,m,n):
	total=m+n+k
	k_frac=float(k)/total
	m_frac=float(m)/total
	n_frac=float(n)/total	
	# dominant - 100% k X k, 0% n X n, 75% m X m, 50% n X m, 100% k X n, 100% k X m
	# opposite: 100% nXn, 25% mXm, 50% nXm
	# m X m: 75%
	#    A    a
	#A  AA   Aa
	#a  Aa   aa
	# m X n: 50%
	#    A    a
	#a  Aa   aa
	#a  Aa   aa
	
	p=(n/total)*((n-1)/(total-1)) + (m/total)*((m-1)/(total-1))*0.25 + (m/total)*(n/(total-1))*0.5 + (n/total)*(m/(total-1))*0.5 
	return 1-p


def main():
	k,m,n=29.0,29.0,18.0
	print get_dominant(k,m,n)


if __name__ == "__main__":
    main()

