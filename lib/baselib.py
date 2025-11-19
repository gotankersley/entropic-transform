def lpad(digits, totalSize):
	itCount = totalSize-len(digits)
	for k in range(itCount):
		digits.insert(0, 0)
	return digits

def b2n(b, n):
	x = 0
	n = n[::-1]
	for i,p in enumerate(n):
		x = x + ((b**i)*p)	  
	return x
   
def n2b(n, b):
	if n == 0:
		return [0]
	digits = []
	while n:
		digits.append(int(n % b))
		n //= b
	return digits[::-1]
	
#https://sochima.me/set-shaping-theory
# compute the information content of a given string
def I(S):
	N = len(S)
	prob = [S.count(c)/N for c in S] # probability of each symbol in the string
	return round(-sum([log2(p) for p in prob]), 2)