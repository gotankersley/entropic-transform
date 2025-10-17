from math import factorial, log2, ceil


class Alpha:
	#Warning: There are several subtle conversions in the alpha <-> perm bijective mapping, and it's easy to get confused
	#Example:
	# - seq = 7247
	# - vals = 2477
	# - alpha = 7724
	# - valById = 0->7,1->2,2->4
	# - idByVal = 7->0,2->1,4->2
	# - repsByVal = 2->1, 4->1,7->2 
	# - repsById = 211 
	# - perm = 0120
	def __init__(self, vals, rot):
		self.rot = rot		
		self.vals = vals[rot:] + vals[:rot] #rotate array around first sequence item	
		self.valById = []
		self.idByVal = {}
		self.repsByVal = {}
		self.repsById = []
		
		id = 0
		prev = -1		
		#Build mapping from vals to permutation id's
		for val in self.vals:
			if val not in self.repsByVal: self.repsByVal[val] = 1
			else: self.repsByVal[val] += 1			
			if val == prev: continue				
			
			self.idByVal[val] = id		
			self.valById.append(val)			
			prev = val			
			id += 1	
			
		for i in range(id):
			val = self.valById[i]
			reps = self.repsByVal[val]
			self.repsById.append(reps)
		
		
	def display(self):
		print ('alpha')
		print ('-rot', self.rot)
		print ('-vals', self.vals)
		print ('-idByVal', self.idByVal)
		print ('-valById', self.valById)
		print ('-repsByVal', self.repsByVal)
		print ('-repsById', self.repsById)
		

def multi_count(n, multiplicities):
	denomProduct = 1
	for m in multiplicities:
		denomProduct *= factorial(m)	
	return factorial(n)//denomProduct
		
def seq_to_multi_perm(seq): 
	vals = seq[:]
	vals.sort() #vals are unordered, but sorted and stored in lexographic order for convenience		
	alpha = Alpha(vals, 0)
	
	perm = []
	for val in seq:
		perm.append(alpha.idByVal[val])
				
	return (perm, alpha) 	
	
	
def multi_perm_to_seq(perm, alpha):
	#Replace permutation id's with actual vals to recreate sequence
	seq = []		
	for id in perm:
		seq.append(alpha.valById[id])
	return seq	
	
		

def multi_rank_to_perm(rank, vals):				
	alpha = Alpha(vals, 0)
	
	#This is the magic method that is at the heart of the compression - unranking a multiset permutation...
	perm = []
	reps = alpha.repsById[:]
	rank += 1 #normalize concepts of zero-based indexing to a 1-based counting by adding 1 to rank
	n = len(alpha.vals)
	
		
	for level in range(n):
		numerator = n-level-1
		numeratorFactorial = factorial(numerator)
		count = 0
		countPrev = 0
		
		for branch,rep in enumerate(reps):
			
			if rep <= 0: continue
			
			#Calculate count of all terminal nodes under current branch
			reps[branch] -= 1 #temporarily decrement as this node is not counted
			#Closed form will be a fraction like 5!/3!1!1! - try to cancel before calculating factorials
			denomProduct = 1
			for r in reps:
				denomProduct *= factorial(r)
			
			count += (numeratorFactorial//denomProduct)		
			reps[branch] += 1 #re-increment to use again
				
			#TODO handle == case?
			if count >= rank:				
				rank -= countPrev				
				id = branch
				perm.append(id)				
				if reps[id] > 0: reps[id] -= 1 #descend down the tree
				break
			else:				
				countPrev = count
			
		
		if count == 0:			
			perm.append(0) #no branches greater than this
			#descend?
	return perm, alpha	
	
		
def multi_perm_to_rank(perm, alpha): #Returns multiset rank
	rank = 0
	
	reps = alpha.repsById[:]
	n = len(alpha.vals)
	
	#Loop over each item in the permutation - each one corresponds to a level in the decision tree
	for level,id in enumerate(perm):
		numerator = n-level-1
		numeratorFactorial = factorial(numerator)
		
		#Loop over all the branches before current branch - find sum of terminal nodes via closed form expression
		for branch in range(id):
			if reps[branch] > 0:
				# From Stack Overflow:
				# "...the number terminal nodes below the current node is equal to the number of permutations using the unused elements in the sequence"
				reps[branch] -= 1 #temporarily decrement as this node is not counted
				#Closed form will be a fraction like 5!/3!1!1! - try to cancel before calculating factorials
				denomProduct = 1
				for r in reps:
					#if r > 0: denomProduct *= r #watch out for zeros = 1												
					denomProduct *= factorial(r)
				#frac = Fraction(numerator, denomProduct) #try to reduce fraction before taking factorial											
				#rank += (factorial(frac.numerator)//factorial(frac.denominator))
				rank += (numeratorFactorial//denomProduct)
				
				reps[branch] += 1 #re-increment to use again

		if reps[id] > 0: reps[id] -= 1 #descend down the tree, decrement node used from available set		
		
	return rank

		


