from math import comb #Binomial coefficient


#NOTE: This ranks/unranks in CO-LEXICAL order, which is different than LEXICAL, but which is easier to generate)
def vals_to_comb_rank(c):
	#CAGES - Algorithm 2.9
	#NOTE: Three important requirements:
	# 1. c must be in decreasing order (i.e., largest values first)
	# 2. All values in c MUST be unique - no duplicates allowed
	# 3. The value zero is not allowed to be in c
	k = len(c)
	r = 0
	for i in range(1, k+1):		
		r = r + comb(c[i-1] - 1, k + 1 - i)
	return r

def vals_to_comb_rank2(c): #Zero index	
	k = len(c)
	r = 0
	for i in range(0, k):	
		r += comb(c[k-1-i], k - i)
	return r
	
	
def comb_rank_to_vals(r, n, k):	
	#CAGES - Algorithm 2.10
	#returns c, the subset of rank, in decreasing order	
	x = n
	c = [0] * k	
	for i in range(1, k+1):  
		while comb(x, k+1-i) > r:    
			x -= 1    
		c[i-1] = x+1
		r = r - comb(x, k + 1 - i)	
	return c
	
	
def comb_rank_to_vals2(r, n, k):#Zero index
	x = n
	c = [0] * k	
	for i in range(0, k):  
		while comb(x, k-i) > r:    
			x -= 1    
		c[k-i-1] = x
		r = r - comb(x, k - i)		
	return c
