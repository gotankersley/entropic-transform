#NOTE: This assumes ZERO indexing on the permutations
#https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.216.9916&rep=rep1&type=pdf
#https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation#C
		
def seq_to_perm(seq):
	vals = sorted(seq)
	return seqvals_to_perm(seq, vals)
	
def seqvals_to_perm(seq, vals): #Perm will have n elements, with unique values from [0 - (n-1)]
	#Assign each value a unique id starting with zero and going to n-1
	vids = {}
	vid = 0
	for v in vals:
		if v not in vids: vids[v] = vid
		vid += 1
	
	perm = []
	dupCounts = {}
	for s in seq:
		if s not in dupCounts:
			perm.append(vids[s])
			dupCounts[s] = 1
		else:
			perm.append(vids[s] + dupCounts[s])
			dupCounts[s] += 1
	
	return perm		
		
		
def permvals_to_seq(perm, vals):
	#Apply the permutation to the values to get a sequence back
	seq = []
	vids = []
	for v in vals:
		vids.append(v)
		
	for p in perm:		
		seq.append(vids[p])
	return seq
	

def myrvold_rank_to_perm(rank, idPerm):	
	#Start with identity so that unranking can shuffle into place		
	n = len(idPerm)
	permRef = idPerm
	unrankRecur(rank, n, permRef) #Output will be placed in permRef
	return permRef
		
		
def perm_to_myrvold_rank(perm): #Returns myrvold rank
	n = len(perm)
	invPerm = [0] * n
	for i in range(n):
		newPos = perm[i]
		invPerm[newPos] = i
		
	permMut = perm[:]
	rank = rankRecur(n, permMut, invPerm)	
	return rank
				
		
		
#Helper functions
def unrankRecur(rank, n, perm):	#Output will be placed in perm
	if n < 1: return #Anchor
	q = rank//n
	r = rank % n
	perm[r], perm[n-1] = perm[n-1], perm[r] #swap
	unrankRecur(q, n-1, perm) #recurse	
	

def rankRecur(n, perm, invPerm): #Returns myrvold rank
	if n < 2: return 0; #Anchor
	
	s = perm[n-1]	
	perm[n-1], perm[invPerm[n-1]] = perm[invPerm[n-1]], perm[n-1] #swap
	invPerm[s], invPerm[n-1] = invPerm[n-1], invPerm[s] #swap
	
	return s + n * rankRecur(n-1, perm, invPerm) #recurse

def getIdentity(n):	
	idPerm = []
	for i in range(n): 
		idPerm.append(i)	
	return idPerm

