from pyscript import window, document
from pyscript.ffi import to_js
from math import floor, ceil, sqrt, gcd, comb, factorial, log2
from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/
from functools import lru_cache
from random import randint
from scipy.special import stirling2
from queue import PriorityQueue

INVALID = -1

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
	

	
#Bijective BWT (BWTS - for BWT "Scottified" by David Scott)
#https://en.wikipedia.org/wiki/Monoid_factorisation#Chen%E2%80%93Fox%E2%80%93Lyndon_theorem
def chen_fox_lyndon_factorization(s: list[int]) -> list[int]:
	"""Monoid factorisation using the Chen–Fox–Lyndon theorem.

	Args:
		s: a list of integers

	Returns:
		a list of integers
	"""
	n = len(s)
	factorization = []
	i = 0
	while i < n:
		j, k = i + 1, i
		while j < n and s[k] <= s[j]:
			if s[k] < s[j]:
				k = i
			else:
				k += 1
			j += 1
		while i <= k:
			factorization.append(s[i:i + j - k])
			i += j - k
	return factorization

from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/
def matrixCompare(a, b):
	itemsA = a['row']
	itemsB = b['row']
	for i in range(len(itemsA)):
		if itemsA[i] > itemsB[i]: return 1
		elif itemsA[i] < itemsB[i]: return -1
		#else continue
	return 0

def digitCompare(a, b):
	 for i in range(len(a)):
		 if a[i] > b[i]: return 1
		 elif a[i] < b[i]: return -1
		 #else continue
	 return 0	
	
def flatten(xss):
	return [x for xs in xss for x in xs]
	
def bwts_transform(n, b, pad): #Bijective Burrows-Wheeler Transform	
	#Note - experimental only.	For serious use, use a suffix array implementation
	seq = n2b(n, b)
	lpad(seq, pad)
	seqLen = len(seq)
	#Example: Banana ->	 ([B], [AN], [AN], [A])
	lyndon = chen_fox_lyndon_factorization(seq)
	#print('Lyndon', lyndon)
	
	matrix = []
	for word in lyndon:	
		wordLen = len(word)
		repWord = (word * ceil(seqLen/wordLen))[:seqLen]
		matrix.append({'row':repWord, 'w':word[-1]})
		
		#Cicular rotations
		for w in range(1, wordLen): #Skip initial rotation
			rotWord = word[w:]+word[:w] #Circular rotatation
			
			#Repeat
			repWord = (rotWord * ceil(seqLen/wordLen))[:seqLen]
			matrix.append({'row':repWord, 'w':rotWord[-1]})	
	
	#Sort - Note we have to sort full matrix with repeated words
	matrix = sorted(matrix, key=cmp_to_key(matrixCompare))	
	
	
	#Get last element in lyndon word - Note, this is NOT just the last column of the matrix
	bwts = []
	for i in range(seqLen):
		bwts.append(matrix[i]['w'])
		#print(matrix[i]['row'])
		
	#dist = b2n(b, bwts)
	return bwts	
	
def bwts_untransform(n, b, pad):	
	#Note - this is an experimental implemention
	#Serious use should use a suffix array implementation instead
	bwts = n2b(n, b)
	lpad(bwts, pad)
	
	seqLen = len(bwts)
	
	#Pivot
	last = []
	for b in bwts:
		last.append([b])		
			
	#Build the rows in reverse (clever trick - see docs for more info on how this actually works)
	rows = [[]]*seqLen	
	for s in range(seqLen): #Sequence iterations
		#print('Iteration:', s, rows)
		for r in range(seqLen): #Row iterations
			rows[r] = last[r] + rows[r]						
		rows = sorted(rows, key=cmp_to_key(digitCompare))
			
		
	seqWords = []
	for r in rows:
		#Check for longest repeated word
		lyndon = chen_fox_lyndon_factorization(r)
		#print('Row', r, 'lyndon', lyndon)
		maxLength = -1
		maxLengthIndex = -1
		w = 0
		for word in lyndon:
			wordLen = len(word)
			if wordLen > maxLength:
				maxLength = wordLen
				maxLengthIndex = w
			w += 1
			
		curWord = lyndon[maxLengthIndex]
		#Check to see if it's a rotation, if so, we can ignore it		
		isValid = True
		for c in range(len(curWord)):
			if curWord[c] != r[c]: #Rotation, so we can ignore
				isValid = False
				break 
				
			
		if isValid:	
			#print('Seq word found', curWord)
			seqWords.append(curWord)		
	seqWords.reverse()
	flatList = flatten(seqWords)
	#dist = b2n(b, flatList)
	return flatList	
	
def seq_to_mtf(seq):
	queue = []
	mtf = []
	largest = max(seq)+1 
	
	for i in range(largest): #Populate queue with index locations
		queue.append(i)
		
	for s in seq:			
		m = queue.index(s)
		mtf.append(m)
		if m != 0: #Update queue
			queue = [queue[m]]+queue[:m]+queue[m+1:] #use array splicing to move-to-front: element-m + elements-before-m + elements-after-m														
		
	return mtf
	
def mtf_to_seq(mtf):
	queue = []
	seq = []
	largest = max(mtf)+1	
	for i in range(largest): #Populate queue with index locations
		queue.append(i)
		
	for m in mtf:		
		seq.append(queue[m])
		if m != 0: #Update queue
			queue = [queue[m]]+queue[:m]+queue[m+1:] #use array splicing to move-to-front: element-m + elements-before-m + elements-after-m								
		
	return seq
## End BWTS section


## Entropic section
#def stirling2(n, k):
#	s = 0
#	for j in range(k+1):			
#		s += (((-1)**(k-j))*(j**n))/(factorial(k-j)*factorial(j))
#	return int(s)

def setPartCompare(a, b):	
	if a[0] > b[0]: return 1
	elif a[0] < b[0]: return -1
	else: return 0	
	
def vals_to_comb_rank2(c): #Zero index	
	k = len(c)
	r = 0
	for i in range(0, k):	
		r += comb(c[k-1-i], k - i)
	return r

def comb_rank_to_vals2(r, n, k):#Zero index
	x = n
	c = [0] * k	
	for i in range(0, k):  
		while comb(x, k-i) > r:    
			x -= 1    
		c[k-i-1] = x
		r = r - comb(x, k - i)		
	return c

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
	
def countSymbolSection(k, n, s): 
	if s < 0 or s > k or s > n: return 0	
	return comb(n, s) * factorial(s) * stirling2(k, s, exact=True)	


def seq_to_set_part(symSeq, symCount):
	setPart = []	
	for i in range(symCount):
		setPart.append([])
			
	for pos in range(len(symSeq)):
		sym = symSeq[pos]
		setPart[sym].append(pos)
		
	#Sort in canonical order
	return sorted(setPart, key=cmp_to_key(setPartCompare))	

@lru_cache(None)
def stirling2MaxLessThan(n, k, m):	
	"""Integer-safe computation of S_{<=m}(n,k) = number of partitions of n labeled items
	into k blocks, each of size <= m."""
	if k == 0:
		return 1 if n == 0 else 0
	if n == 0:
		return 1 if k == 0 else 0
	# DP table: S[i][j] = S_{<=m}(i,j)
	S = [[0]*(k+1) for _ in range(n+1)]
	S[0][0] = 1
	for ni in range(1, n+1):
		for ki in range(1, min(k, ni)+1):
			total = 0
			for j in range(1, min(m, ni)+1):
				total += comb(ni-1, j-1) * S[ni-j][ki-1]
			S[ni][ki] = total	
	return S[n][k]
		
	
#This calculates the Stirling2 number of set-parts with n elements, k parts, and having a max part size that is exactly m
def stirling2Max(n, k, m):
	return stirling2MaxLessThan(n, k, m) - stirling2MaxLessThan(n, k, m-1)
	
#This calculates the Stirling2 number of set-parts with n elements, k parts, and having a max part size that is exactly m, and has an initial part size of r
def stirling2MaxR(n, k, m, r):
	cmb = comb(n-1,r-1)
	if r == m:
		#TODO - For some reason this is the total minus all the other ones
		total = stirling2Max(n, k, m)
		previousSum = 0
		for i in range(1, r):
			previousSum += stirling2MaxR(n, k, m, i)
		return total-previousSum
	else: return (cmb*stirling2MaxLessThan(n-r,k-1, m)) - (cmb*stirling2MaxLessThan(n-r,k-1, m-1))

#Rank an RGF given as a list/array of integers (1-based labels).
#Returns 0-based lexicographic rank among all RGFs of the same (n,k).
def stirling_rank(rgf_list):

	rgf = list(rgf_list)
	n = len(rgf)
	if n == 0:
		return 0
	if rgf[0] != 1:
		raise ValueError("RGF must start with 1 (use 1-based labels).")
	k = max(rgf)
	max_so_far = 1
	for i in range(1, n):
		allowed = max_so_far + 1
		if not (1 <= rgf[i] <= allowed):
			raise ValueError(f"Invalid RGF at position {i}: {rgf[i]} (allowed 1..{allowed})")
		if rgf[i] > max_so_far:
			max_so_far = rgf[i]

	@lru_cache(None)
	def dp(pos, used, maxval): #dp - Dynamic Programming
		if pos == n:
			return 1 if used == k else 0
		total = 0
		for val in range(1, maxval+2):
			new_used = used if val <= maxval else used+1
			if new_used <= k:
				total += dp(pos+1, new_used, max(maxval, val))
		return total

	rank = 0
	used = 1
	maxval = 1
	for pos in range(1, n):
		current = rgf[pos]
		for smaller in range(1, current):
			new_used = used if smaller <= maxval else used+1
			if new_used <= k:
				rank += dp(pos+1, new_used, max(maxval, smaller))
		if current > maxval:
			used += 1
			maxval = current
	return rank
	

#Unrank to an RGF list (1-based labels) of length n with k blocks.
#rank is 0-based.
def stirling_unrank(n, k, rank):
	from functools import lru_cache
	@lru_cache(None)
	def dp(pos, used, maxval): #dp - Dynamic Programming
		if pos == n:
			return 1 if used == k else 0
		total = 0
		for val in range(1, maxval+2):
			new_used = used if val <= maxval else used+1
			if new_used <= k:
				total += dp(pos+1, new_used, max(maxval, val))
		return total

	total_rgfs = dp(1,1,1)
	if rank < 0 or rank >= total_rgfs:
		raise ValueError(f"rank out of range: 0..{total_rgfs-1}")

	rgf = [1]
	used, maxval = 1, 1
	for pos in range(1, n):
		for val in range(1, maxval+2):
			new_used = used if val <= maxval else used+1
			if new_used > k:
				continue
			cnt = dp(pos+1, new_used, max(maxval, val))
			if rank < cnt:
				rgf.append(val)
				if val > maxval:
					used += 1
					maxval = val
				break
			rank -= cnt
	return rgf
	
def rgf_to_set_part(rgf):
	"""
	Convert an RGF list (1-based labels) into a set-partition
	represented as a list of blocks (each block is a list of element indices).
	
	Example:
		[1,1,2,2] -> [[0,1],[2,3]]
	"""
	blocks = {}
	for idx, label in enumerate(rgf):
		blocks.setdefault(label, []).append(idx)
	# return blocks in label order
	return [blocks[k] for k in sorted(blocks)]


def set_part_to_rgf(setPart):#, n=None):
	"""
	Convert a set-partition (list of blocks) into an RGF list (1-based labels).
	Each block is a list of element indices.
	
	Example:
		[[0,1],[2,3]] -> [1,1,2,2]
	"""
	#if n is None:
	n = sum(len(block) for block in setPart)
	rgf = [0] * n
	for label, block in enumerate(setPart, start=1):
		for idx in block:
			rgf[idx] = label
	return rgf
	
def getSizeOfLargestPart(setPart):
	largestPart = -1
	for p in setPart:
		largestPart = max(largestPart, len(p))
	return largestPart

def getValsFromCombVals(combVals): 
	combVals.reverse()
	for i in range(len(combVals)):
		combVals[i] -= 1
	return combVals

def getValToSymMap(vals):
	val_to_sym = {}	
	
	sym = 0
	for val in vals:
		val_to_sym[val] = sym				
		sym += 1
		
	return val_to_sym 
	
def getSymToValMap(vals):	
	sym_to_val = []
	
	sym = 0
	for val in vals:		
		sym_to_val.append(val)
		sym += 1
		
	return sym_to_val
	


def getSymSeq(valSeq, val_to_sym):	
	symSeq = []
	for val in valSeq:		
		symSeq.append(val_to_sym[val])
	return symSeq


def getSymPerm(symSeq, setPart):
	symPerm = []
	for part in setPart:
		#Look at initial to get pos in seq
		initialPos = part[0]
		symPerm.append(symSeq[initialPos])
	return symPerm

def getCombVals(vals):
	#COLEX rank must be in descending order
	combVals = []
	for i in range(len(vals)-1, -1, -1):
		combVals.append(vals[i]+1)
	return combVals
	
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
	
		
	seq = near_entropic_unrank(n, b, pad)
	dist = b2n(b, seq)
	return dist
	
def entropic_rank(valSeq, totalSymbolsAvail):	
	seqLen = len(valSeq)		
	
	vals = sorted(list(set(valSeq)))
	symCount = len(vals)		
	
	val_to_sym = getValToSymMap(vals)
	
	symSeq = getSymSeq(valSeq, val_to_sym) 		
	
	rank = 0

	
	#  1. Add symbol sections
	for i in range(1, symCount):
		rank += countSymbolSection(seqLen, totalSymbolsAvail, i)
	#Calculate section sizes	
	totalComb = comb(totalSymbolsAvail, symCount)
	totalSymPerm = factorial(symCount)
		
	stirSectionSize = totalComb * totalSymPerm
	combSectionSize = totalSymPerm	

	
	#NOTE: Ranking the set-partition is the real heart of this algorithm.  
	#As the CAGES book says on ranking Set Partitions: "It is both convenient and natural to use an RGF for ranking".
	#Unfortunately, this method doesn't use an RGF, (and, thus, is neither convenient or natural).  
	#However, the extra complexity is required in order to have it rank in fully entropic order, as the RGF order is not entropic.
	setPart = seq_to_set_part(symSeq, symCount)	
	n = seqLen	
	stirRank = 0
	prevLargestPartSize = (n-symCount) + 1
	usedElements = list(range(n))
				
	for p in range(symCount-1): #Minus-1, because there is only a single way to do the last part that uses all  the remaining elements
		curSetPart = setPart[p:]		
		k = len(curSetPart)

			
		
		#  2. By Set Partition largest part section					
		m = getSizeOfLargestPart(curSetPart)
		prevStir = stirRank		
		for i in range(prevLargestPartSize, m, -1):
			largestPartSectionSize = stirling2Max(n, k, i)			
			stirRank += largestPartSectionSize
		prevLargestPartSize = m
		#Note: At this point in the unranking we have just found the value of m (largest part size)
					
		
		#  3. By Set Partition initial part size section			
		r = len(curSetPart[0])
		for i in range(m, r, -1):			
			initialPartSectionSize = stirling2MaxR(n, k, m, i)			
			stirRank += initialPartSectionSize #Some of these will be zero, but they'll just be ignored
		#Note: At this point in the unranking we know the length of the initial part of the set partition r		
		initialPartSectionSize = stirling2MaxR(n, k, m, r)										
		
		#  4. By Set Partition initial element combination - (The first element is always fixed as the lowest)
		elementCombs = comb(n-1, r-1)
		elementSectionSize = initialPartSectionSize // elementCombs
			
		elementIds = []		
		usedElements.pop(0)
		if r == 1: pass
		else:
			usedElementsCopy = usedElements[:]
			for i in range(1, r):
				v = curSetPart[0][i]
				id = usedElementsCopy.index(v)
				id2 = usedElements.index(v)
				elementIds.append(id)
				usedElements.pop(id2)			
			elementRank = vals_to_comb_rank2(elementIds) 	
			stirRank += (elementSectionSize * elementRank)	
		
		#Iterate
		n -= r
		
		
	rank += (stirRank * stirSectionSize)
		
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy, 
	#         so the rest is just one possible scheme to identify the sequence	
	#  -------------
	
	
	#  5. Add the combination rank 	
	combRank = vals_to_comb_rank2(vals)
	rank += (combSectionSize * combRank)	
	
	
	#  6. Add the Sym/Set-Part Perm Rank (Myrvold)	
	symPerm = getSymPerm(symSeq, setPart)
	symRank = perm_to_myrvold_rank(symPerm)	
	rank += symRank


	return rank
	
	
def entropic_unrank(rank, totalSymbolsAvail, seqLen):		
			
		
	#  1. Get symbol section	
	count = 0
	prevCount = count
	for i in range(1, totalSymbolsAvail+1):
		count += countSymbolSection(seqLen, totalSymbolsAvail, i)
		if count > rank:
			symCount = i			
			rank -= prevCount
			break
		else: prevCount = count
	try: symCount
	except: symCount = totalSymbolsAvail
		
	#Calculate section sizes	
	totalComb = comb(totalSymbolsAvail, symCount)
	totalSymPerm = factorial(symCount)
			
	stirSectionSize = totalComb * totalSymPerm
	combSectionSize = totalSymPerm	
	
	# Unrank the Set-Partition one part at a time
	stirRank = rank // stirSectionSize	
	setPart = [] 
	n = seqLen
	usedElements = list(range(0,n))
	prevLargestPartSize = (n-symCount)+1

	for p in range(symCount-1): #Minus-1, because there is only a single way to do the last part that uses all  the remaining elements		
		k = symCount - p
		
		#  2. By Set Partition largest part section	
		count = 0		
		prevCount = 0		
		#for m in range((n-k) + 1, 0, -1): 			
		for m in range(prevLargestPartSize, 0, -1): 			
			largestPartSectionSize = stirling2Max(n, k, m)			
			count += largestPartSectionSize			
			if count > stirRank: 
				stirRank -= prevCount
				break
			else: prevCount = count		
		try: m
		except: m = (n-k)+1		
		prevLargestPartSize = m

		
		#  3. By Set Partition initial part size section	
		count = 0
		prevCount = 0		
		for r in range(m, 0, -1):			
			initialPartSectionSize = stirling2MaxR(n, k, m, r) #Some of these will be zero, but they'll just be ignored
			count += initialPartSectionSize						
			if count > stirRank: 				
				stirRank -= prevCount
				#initialPartSectionSize = count
				break
			else: prevCount = count
		initialPartSectionSize = stirling2MaxR(n, k, m, r)	
		if initialPartSectionSize == 0: initialPartSectionSize = 1

		
		
		#  4. By Set Partition initial element combination - (The first element is always fixed as the lowest)		
		elementCombs = comb(n-1, r-1)		
		elementSectionSize = initialPartSectionSize // elementCombs		
		elementRank = stirRank // elementSectionSize
		elementIds = comb_rank_to_vals2(elementRank, n-1, r-1) 		
		elementsVals = [usedElements.pop(0)]	
		if r == 1: pass
		else:
			usedElementsCopy = usedElements[:]
			for id in elementIds:					
				elementVal = usedElementsCopy[id]			
				id2 = usedElements.index(elementVal)
				usedElements.pop(id2)			
				elementsVals.append(elementVal)										

		setPart.append(elementsVals)
		stirRank -= (elementSectionSize * elementRank)
		
		#Iterate				
		n -= r
			
	#Use all remaining elements for the last set	
	elements = []
	for i in range(len(usedElements)):			
			element = usedElements.pop(0)
			elements.append(element)
	setPart.append(elements)	

	
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy	
	#  -------------
		
	
	##  5. Get the values from the combination rank of symbols	
	combRank = (rank % stirSectionSize) // combSectionSize
	vals = comb_rank_to_vals2(combRank, totalSymbolsAvail, symCount) 	
	
	
	
	#  6. Get the Sym/Set-Part Perm 
	symRank = rank % totalSymPerm
	idPerm = getIdentity(symCount)
	symPerm = myrvold_rank_to_perm(symRank, idPerm)		
	
	
	# Apply to recreate seq
	sym_to_val = getSymToValMap(vals)
	seq = [0] * seqLen
	sym = 0
	for part in setPart:
		val = sym_to_val[symPerm[sym]]
		for pos in part:
			seq[pos] = val
		sym += 1
				
	return seq
	
def near_entropic_rank(valSeq, totalSymbolsAvail):	
	seqLen = len(valSeq)		
	
	vals = sorted(list(set(valSeq)))
	symCount = len(vals)		
	
	val_to_sym = getValToSymMap(vals)
	
	symSeq = getSymSeq(valSeq, val_to_sym) 		
	
	rank = 0
	
	#print('Sym Seq', symSeq)
	#print('Symbols:', symCount)
	
	
	#  1. Add symbol sections
	for i in range(1, symCount):
		rank += countSymbolSection(seqLen, totalSymbolsAvail, i)
	#print('rank after symbol section', rank)
	
	#Calculate section sizes	
	totalComb = comb(totalSymbolsAvail, symCount)
	totalSymPerm = factorial(symCount)
		
	stirSectionSize = totalComb * totalSymPerm
	combSectionSize = totalSymPerm
				
	#  2. Add Set Partition / Stirling2 rank 
	setPart = seq_to_set_part(symSeq, symCount)
	rgf = set_part_to_rgf(setPart)		
	stirRank = stirling_rank(rgf) 
	rank += (stirRank * stirSectionSize)
	#print ('Stir Section size', stirSectionSize)
	#print ('Set Part', setPart)
	#print ('RGF', rgf)
	#print ('Stir Rank', stirRank)
	
		
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy, 
	#         so the rest is just one possible scheme to identify the sequence	
	#  -------------
	
	
	#  3. Add the combination rank 
	combVals = getCombVals(vals)		
	combRank = vals_to_comb_rank(combVals) #COLEX order	
	rank += (combSectionSize * combRank)	
	#print ('Comb Section size: ', combSectionSize)
	#print ('Comb Vals: ', combVals)
	#print('Comb Rank: ', combRank)
	
	
	#  4. Add the Sym/Set-Part Perm Rank (Myrvold)	
	symPerm = getSymPerm(symSeq, setPart)
	symRank = perm_to_myrvold_rank(symPerm)	
	rank += symRank
	#print ('Sym Perm: ', symPerm)	
	#print ('Sym Rank: ', symRank)
		
	
	#print('Rank:', rank)
	return rank
	
	
def near_entropic_unrank(rank, totalSymbolsAvail, seqLen):		
			
	
	#  1. Get symbol section
	#print('Starting Rank', rank)
	count = 0
	prevCount = count
	for i in range(1, totalSymbolsAvail+1):
		count += countSymbolSection(seqLen, totalSymbolsAvail, i)
		if count > rank:
			symCount = i			
			rank -= prevCount
			break
		else: prevCount = count
	try: symCount
	except: symCount = totalSymbolsAvail
	
	#print ('Symbol Section', symCount)
	#print ('Rank after symbol section:', rank)
		
	#Calculate section sizes	
	totalComb = comb(totalSymbolsAvail, symCount)
	totalSymPerm = factorial(symCount)
		
	stirSectionSize = totalComb * totalSymPerm
	combSectionSize = totalSymPerm
	
	#  2. Get the Set Partition from Stirling2 rank
	stirRank = rank // stirSectionSize
	rgf = stirling_unrank(seqLen, symCount, stirRank)
	setPart = rgf_to_set_part(rgf)		
	#print ('Stir Rank', stirRank)	
	#print ('RGF', rgf)
	#print ('Set Part', setPart)

	
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy	
	#  -------------
		
	
	#  3. Get the values from the combination rank of symbols	
	combRank = (rank % stirSectionSize) // combSectionSize
	combVals = comb_rank_to_vals(combRank, totalSymbolsAvail, symCount) #COLEX order
	combVals.reverse()
	for i in range(len(combVals)):
		combVals[i] -= 1
	vals = combVals
	#print ('Comb Section size: ', combSectionSize)
	#print ('Comb Rank: ', combRank)
	#print ('Comb Vals: ', combVals)
	
	
	#  4. Get the Sym/Set-Part Perm 
	symRank = rank % totalSymPerm
	idPerm = getIdentity(symCount)
	symPerm = myrvold_rank_to_perm(symRank, idPerm)	
	#print ('Sym Perm: ', symPerm)	
	#print ('Sym Rank: ', symRank)
	
	
	# Apply to recreate seq
	sym_to_val = getSymToValMap(vals)
	seq = [0] * seqLen
	sym = 0
	for part in setPart:
		val = sym_to_val[symPerm[sym]]
		for pos in part:
			seq[pos] = val
		sym += 1
			
	#print('Seq:', seq)
	return seq
## End Entropic Section
	
	
def I(S):
	N = len(S)
	prob = [S.count(c)/N for c in S] # probability of each symbol in the string
	return round(-sum([log2(p) for p in prob]), 2)
	
def calculate(e):
	global menu
	menu = window.menu
		
	seqLen = int(document.querySelector('#seqLen').value)
	symCount = int(document.querySelector('#symCount').value)
	
	#Generate Random
	rndSeq = []
	for i in range(seqLen):
		rndSeq.append(randint(0, symCount-1))
			
	rndRank = b2n(symCount, rndSeq)
		
		
	#Output
	data = {		
		'rndSeq':rndSeq, 	  
		'rndRank':rndRank, 
		'rndEntropy':I(rndSeq), 
	}
	
	
	if menu.nearerEntropic:
		entRank = entropic_rank(rndSeq, symCount)
		entShapedSeq = entropic_unrank(entRank, symCount, seqLen + menu.kSetShaping)		
		data['entShapedSeq'] = entShapedSeq
		data['entRank'] = entRank
		data['entShapedEntropy'] = I(entShapedSeq)
		
		entSeqBaseA = entropic_unrank(rndRank, symCount, seqLen + menu.kSetShaping)
		data['entSeqBaseA'] = entSeqBaseA	
		data['entEntropyBaseA'] = I(entSeqBaseA)
		
	if menu.nearEntropic:
		nearEntRank = near_entropic_rank(rndSeq, symCount)
		nearEntShapedSeq = near_entropic_unrank(nearEntRank, symCount, seqLen + menu.kSetShaping)		
		data['nearEntShapedSeq'] = nearEntShapedSeq
		data['nearEntRank'] = nearEntRank
		data['nearEntShapedEntropy'] = I(nearEntShapedSeq)
		
		nearSeqBaseA = near_entropic_unrank(rndRank, symCount, seqLen + menu.kSetShaping)
		data['nearSeqBaseA'] = nearSeqBaseA		
		data['nearEntropyBaseA'] = I(nearSeqBaseA)
		
	if menu.bwts:
		bwtsSeq = bwts_transform(rndRank, symCount, seqLen)		
		bwtsMtfSeq = seq_to_mtf(bwtsSeq)
		data['bwtsMtfSeq'] = bwtsMtfSeq
		data['bwtsMtfEntropy'] = I(bwtsMtfSeq)
		data['bwtsMtfRank'] = b2n(symCount, bwtsMtfSeq)
		
	
	window.showResults(to_js(data))