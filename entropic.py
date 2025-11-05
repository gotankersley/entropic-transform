from scipy.special import stirling2
from math import comb, factorial

from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/

from lib.myrvold import *
from lib.combinations import *
from lib.stirling import *

DEBUG = False


	
#ABOUT:  This is a function to rank a sequence in sorted entropic order, rather than the regular place-value base system.
# This means that sequences with the less entropy, (i.e. fewer distinct symbols), like 'AAA', 'BBB', are ranked first,
# and those with the more entropy, (i.e. more differing symbols), like 'ABCD', are ranked last.
#
# Sequences are uniquely identified in the following way:
#  1. By symbol section, (i.e. all the sequences with only 1 symbol are grouped together, and so on)
#  2. By Set Partition largest part section			--|
#  3. By Set Partition initial part section			  +-- 2,3,4 Are in lieu of an RGF
#  4. By Set Partition initial element combination	--|
#  -------------	
#	Note: Everything after this has the exact same amount of entropy, 
#         so the rest is just one possible scheme to identify the sequence
#  -------------
#  5. By the combination rank of the symbols used from the set of total available symbols
#  6. By the permutation rank of the mapping the symbols to the set partition
def entropic_rank(valSeq, totalSymbolsAvail):	
	seqLen = len(valSeq)		
	
	vals = sorted(list(set(valSeq)))
	symCount = len(vals)		
	
	val_to_sym = getValToSymMap(vals)
	
	symSeq = getSymSeq(valSeq, val_to_sym) 		
	
	rank = 0
	
	if DEBUG:
		print('Starting Ranking:')
		print('\t', '- Val Seq: ', valSeq)
		print('\t', '- Sym Seq', symSeq)
		print('\t', '- Sym Count:', symCount)
		print()
	
	#  1. Add symbol sections
	for i in range(1, symCount):
		rank += countSymbolSection(seqLen, totalSymbolsAvail, i)
	#Calculate section sizes	
	totalComb = comb(totalSymbolsAvail, symCount)
	totalSymPerm = factorial(symCount)
		
	stirSectionSize = totalComb * totalSymPerm
	combSectionSize = totalSymPerm	
	if DEBUG: 
		print('1. Symbol sections')
		print('\t', '- Rank after symbol section', rank)
		print('\t', '- TotalComb', totalComb)
		print('\t', '- TotalSymPerm', totalSymPerm)
		print('\t', '- StirSectionSize:', stirSectionSize)
		print('\t', '- CombSectionSize:', combSectionSize)		
		print()
	
	#NOTE: Ranking the set-partition is the real heart of this algorithm.  
	#As the CAGES book says on ranking Set Partitions: "It is both convenient and natural to use an RGF for ranking".
	#Unfortunately, this method doesn't use an RGF, (and, thus, is neither convenient or natural).  
	#However, the extra complexity is required in order to have it rank in fully entropic order, as the RGF order is not entropic.
	setPart = seq_to_set_part(symSeq, symCount)	
	n = seqLen	
	stirRank = 0
	prevLargestPartSize = (n-symCount) + 1
	usedElements = list(range(n))
	if DEBUG:
		print('2. Set Part Rank:', setPart)		
		print('\t', '- usedElements:', usedElements)				
	for p in range(symCount-1): #Minus-1, because there is only a single way to do the last part that uses all  the remaining elements
		curSetPart = setPart[p:]		
		k = len(curSetPart)
		if DEBUG:
			print('SET-PART RANK ITERATION: ', p)
			print('--', 'Cur Set Part', curSetPart)
			print('--', 'k=', k)
			print()
			
		
		#  2. By Set Partition largest part section					
		m = getSizeOfLargestPart(curSetPart)
		prevStir = stirRank		
		for i in range(prevLargestPartSize, m, -1):
			largestPartSectionSize = stirling2Max(n, k, i)			
			stirRank += largestPartSectionSize
		prevLargestPartSize = m
		#Note: At this point in the unranking we have just found the value of m (largest part size)
		if DEBUG:			
			largestPartSectionSize = stirling2Max(n, k, m)			
			print('--', 'Stir Rank before Largest part section size: ', prevStir)	
			print('--', 'Stir Rank after Largest part section size: ', stirRank)	
			print('--', 'Section size: ', largestPartSectionSize)
			print('--', 'm=', m)
			print()
				
		
		#  3. By Set Partition initial part size section			
		r = len(curSetPart[0])
		for i in range(m, r, -1):			
			initialPartSectionSize = stirling2MaxR(n, k, m, i)			
			stirRank += initialPartSectionSize #Some of these will be zero, but they'll just be ignored
		#Note: At this point in the unranking we know the length of the initial part of the set partition r		
		initialPartSectionSize = stirling2MaxR(n, k, m, r)				
		if DEBUG:			
			print('--', 'Stir Rank after initial part section size: ', stirRank)
			print('--', 'Section size: ', initialPartSectionSize)
			print('--', 'r=', r)
			print()
						
		
		#  4. By Set Partition initial element combination - (The first element is always fixed as the lowest)
		elementCombs = comb(n-1, r-1)
		elementSectionSize = initialPartSectionSize // elementCombs
			
		elementIds = []		
		usedElements.pop(0)
		if r == 1:
			if DEBUG:
				print('--','Part has only single element')
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
			if DEBUG:
				print('--', 'Stir Rank after element section', stirRank)
				print('--', 'ElementSectionSize', elementSectionSize)
				print('--', 'ElementCombs', elementCombs)
				print('--', 'Element Ids', elementIds)
				print('--', 'ElementRank', elementRank)
				print('--', 'Used Elements Before', usedElementsCopy)
				print('--', 'Used Elements After', usedElements)			
				print()
		
		#Iterate
		if DEBUG:
			print('Iterating:', 'New n=', n-r, '\r\n')
		n -= r
		
		
	rank += (stirRank * stirSectionSize)
	if DEBUG:
		print('Final Set Part Rank: ', stirRank)
		print('Rank after adding Set Part Rank: ', rank)
		
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy, 
	#         so the rest is just one possible scheme to identify the sequence	
	#  -------------
	
	
	#  5. Add the combination rank 	
	combRank = vals_to_comb_rank2(vals)
	rank += (combSectionSize * combRank)	
	if DEBUG:
		print('5. Combination rank')
		print('\t', '- Rank after combination section', rank)
		print('\t', '- Comb Vals: ', vals)
		print('\t', '- Comb Rank: ', combRank)
		print('\t', '- Comb Section size: ', combSectionSize)
		print()
	
	
	#  6. Add the Sym/Set-Part Perm Rank (Myrvold)	
	symPerm = getSymPerm(symSeq, setPart)
	symRank = perm_to_myrvold_rank(symPerm)	
	rank += symRank
	if DEBUG:
		print('6. Sym rank')
		print('\t', '- Rank after sym section', rank)
		print ('\t', '- Sym Perm: ', symPerm)	
		print ('\t', '- Sym Rank: ', symRank)
		print()
		
	
	if DEBUG: print('Final Rank:', rank)	

	return rank
	
	
def entropic_unrank(rank, totalSymbolsAvail, seqLen):		
			
	if DEBUG:
		print('Starting Unranking:')
		print('\t', '- Rank: ', rank)
		print('\t', '- Total Symbols Avail', totalSymbolsAvail)
		print('\t', '- Seq Len:', seqLen)
		print()
		
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
	if DEBUG: 
		print('1. Symbol sections')
		print('\t', '- Sym Count: ', symCount)
		print('\t', '- Rank after symbol section', rank, prevCount)
		print('\t', '- TotalComb', totalComb)
		print('\t', '- TotalSymPerm', totalSymPerm)
		print('\t', '- StirSectionSize:', stirSectionSize)
		print('\t', '- CombSectionSize:', combSectionSize)		
		print()
	
	# Unrank the Set-Partition one part at a time
	stirRank = rank // stirSectionSize	
	setPart = [] 
	n = seqLen
	usedElements = list(range(0,n))
	prevLargestPartSize = (n-symCount)+1
	if DEBUG:
		print('2. Set Part Unranking:', setPart)			
		print('\t', '- StirRank:', stirRank)	
		print('\t', '- usedElements:', usedElements)
	for p in range(symCount-1): #Minus-1, because there is only a single way to do the last part that uses all  the remaining elements		
		k = symCount - p
		if DEBUG:
			print('SET-PART ITERATION: ', p)
			print('--', 'Set Part', setPart)
			print('--', 'k=', k)
			print()		
		
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
		if DEBUG:
			largestPartSectionSize = stirling2Max(n, k, m)				
			print('--', 'Stir Rank before Largest part section: ', stirRank+prevCount)
			print('--', 'Stir Rank after Largest part section: ', stirRank)			
			print('--', 'Section size: ', largestPartSectionSize)
			print('--', 'Size preceeding: ', prevCount)
			print('--', 'Unranked m: ', m)
			print()
		
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
		if DEBUG:			
			print('--', 'Stir Rank before initial part section: ', stirRank+prevCount)
			print('--', 'Stir Rank after initial part section size: ', stirRank)
			print('--', 'Section size: ', initialPartSectionSize)
			print('--', 'Unranked r: ', r)
			print()
		
		
		#  4. By Set Partition initial element combination - (The first element is always fixed as the lowest)		
		elementCombs = comb(n-1, r-1)		
		elementSectionSize = initialPartSectionSize // elementCombs		
		elementRank = stirRank // elementSectionSize
		elementIds = comb_rank_to_vals2(elementRank, n-1, r-1) 		
		elementsVals = [usedElements.pop(0)]	
		if r == 1:
			if DEBUG:
				print('--', 'Part has only single element: ', elementsVals)				
		else:
			usedElementsCopy = usedElements[:]
			for id in elementIds:					
				elementVal = usedElementsCopy[id]			
				id2 = usedElements.index(elementVal)
				usedElements.pop(id2)			
				elementsVals.append(elementVal)										
			if DEBUG:
				print('--', 'Stir Rank after element section', stirRank)
				print('--', 'ElementSectionSize', elementSectionSize)
				print('--', 'ElementCombs', elementCombs)
				print('--', 'Element Ids', elementIds)
				print('--', 'Element Vals', elementsVals)
				print('--', 'ElementRank', elementRank)
				print('--', 'Used Elements Before', usedElementsCopy)
				print('--', 'Used Elements After', usedElements)			
				print()
		setPart.append(elementsVals)
		stirRank -= (elementSectionSize * elementRank)
		
		#Iterate		
		if DEBUG:
			print('Iterating:', 'New n=', n-r, '\r\n')
		n -= r
			
	#Use all remaining elements for the last set	
	elements = []
	for i in range(len(usedElements)):			
			element = usedElements.pop(0)
			elements.append(element)
	setPart.append(elements)
	if DEBUG:
		print('Unranked Set Part: ', setPart)	

	
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy	
	#  -------------
		
	
	##  5. Get the values from the combination rank of symbols	
	combRank = (rank % stirSectionSize) // combSectionSize
	vals = comb_rank_to_vals2(combRank, totalSymbolsAvail, symCount) 	
	if DEBUG:
		print('5. Combination rank')
		print('\t', '- Comb Rank: ', combRank)
		print('\t', '- Comb Vals: ', vals)		
		print('\t', '- Comb Section size: ', combSectionSize)
		print()	
	
	
	#  6. Get the Sym/Set-Part Perm 
	symRank = rank % totalSymPerm
	idPerm = getIdentity(symCount)
	symPerm = myrvold_rank_to_perm(symRank, idPerm)	
	if DEBUG:
		print('6. Sym rank')		
		print ('\t', '- Sym Rank: ', symRank)
		print ('\t', '- Sym Perm: ', symPerm)	
		print()	
	
	
	# Apply to recreate seq
	sym_to_val = getSymToValMap(vals)
	seq = [0] * seqLen
	sym = 0
	for part in setPart:
		val = sym_to_val[symPerm[sym]]
		for pos in part:
			seq[pos] = val
		sym += 1
			
	if DEBUG: print('Final Seq:', seq)
	return seq
	

#Combinatorial function to count the ways to rank:
# - Filling up a sequence of k length (e.g. k = 3)
# - Count of set of n symbols available (e.g. n = {A,B,C} = 3 )
# - Using subset s of symbols (e.g. s = 2)
def countSymbolSection(k, n, s): 
	if s < 0 or s > k or s > n: return 0	
	return comb(n, s) * factorial(s) * stirling2(k, s, exact=True)
	
		
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