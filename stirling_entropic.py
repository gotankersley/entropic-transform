from scipy.special import stirling2
from math import comb, factorial

from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/

from lib.stirling import *
from lib.combinations import *

DEBUG = False


	
#ABOUT
# Sequences are uniquely identified in the following way:
#  1. By Set Partition largest part section			--|
#  2. By Set Partition initial part section			  +-- 2,3,4 Are in lieu of an RGF
#  3. By Set Partition initial element combination	--|
def stirling_entropic_rank(setPart, seqLen, symCount):	
	
	#NOTE:
	#As the CAGES book says on ranking Set Partitions: "It is both convenient and natural to use an RGF for ranking".
	#Unfortunately, this method doesn't use an RGF, (and, thus, is neither convenient or natural).  
	#However, the extra complexity is required in order to have it rank in fully entropic order, as the RGF order is not entropic.	
	n = seqLen	
	stirRank = 0
	prevLargestPartSize = (n-symCount) + 1
	usedElements = list(range(n))
	if DEBUG:
		print('1. Set Part Rank:', setPart)		
		print('\t', '- usedElements:', usedElements)	
		
	for p in range(symCount-1): #Minus-1, because there is only a single way to do the last part that uses all the remaining elements
		curSetPart = setPart[p:]		
		k = len(curSetPart)
		if DEBUG:
			print('SET-PART RANK ITERATION: ', p)
			print('--', 'Cur Set Part', curSetPart)
			print('--', 'k=', k)
			print()
			
		
		#  1. By Set Partition largest part section					
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
				
		
		#  2. By Set Partition initial part size section			
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
						
		
		#  3. By Set Partition initial element combination - (The first element is always fixed as the lowest)
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
		
		
	return stirRank
	
	
def stirling_entropic_unrank(stirRank, seqLen, symCount):		
			
	if DEBUG:
		print('Starting Unranking:')
		print('\t', '- Stir Rank: ', stirRank)		
		print('\t', '- Seq Len:', seqLen)
		print()
	
	
	
	# Unrank the Set-Partition one part at a time	
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
		
		#  1. By Set Partition largest part section	
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
		
		#  2. By Set Partition initial part size section	
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
		
		
		#  3. By Set Partition initial element combination - (The first element is always fixed as the lowest)		
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

	return setPart
	

	
		
