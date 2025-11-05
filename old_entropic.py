#***********************
#NOTE:  DEPRECATED!!
#		While this is technically correct, it unfortunately slows down very quickly for increasingly large amounts of symbols, (but not sequence length)
#		For large amounts of symbols, try using the near-entropic rank instead.
#***********************
from scipy.special import stirling2
from math import comb, factorial
from collections import Counter

from lib.combinations import *
from lib.int_partitions import *
from lib.multiset import *


	
#ABOUT:  This is a function to rank a sequence in sorted entropic order, rather than the regular place-value base system.
# This means that sequences with the least entropy, (i.e. fewest distinct symbols), like 'AAA', 'BBB', are ranked first,
# and those with the most entropy, (i.e. most differing symbols), like 'ABCD', are ranked last.
# Now, this is rather tricky, because the most obvious way to rank, 
# (using the recursive suffix repeat counts - e.g. BCAB = [[2,1,1],[1,1,1],[1,1][1]]), requires some count function
# that I wasn't able to find. But, even if it could be found, this isn't enough to uniquely identify a sequence,
# there is still some additional permutations within that, which would still need to be identified.
# Because of that, this is less obvious method that uniquely identifies sequences in the following way:
#  1. By symbol section, (i.e. all the sequences with only 1 symbol are grouped together, and so on)
#  2. By repeat pattern, this is the count that each symbol appears, ordered by size.  (e.g. BCAB = [2,1,1], or BAAA = [3,1])
#  -------------
#   Note: Everything after this has the exact same amount of entropy, 
#         so the rest is just one possible scheme to identify a sequence,
#		  (and that order is similar to, but subtly different than lexical)
#  -------------
#  3. By the combination rank of the count of symbols available, and the size of the repeat patterns
#  4. By the multi-set permutation rank of the repeats (mspReps)
#  5. By multi-set permutation rank of the sequence (mspSeq)
#
#  So yes, 5 different, crazy, individual levels required to uniquely identify, (and rank) a sequence in entropic order!
#  [sym][reps][comb][mspReps][mspSeq]
def entropic_rank(seq, totalSymbolsAvail):		
	k = len(seq)	
	n = totalSymbolsAvail
	
	symbolsUsed = set(seq)
	s = len(symbolsUsed)
	reps = getReps(seq)
	rank = 0
	
	#print('Seq', seq)
	#print('Symbols:', symbolsUsed, s)
	#print('Reps:', reps)
	
	#  1. Add symbol sections
	for i in range(1, s):
		rank += countSymbolSection(k, n, i)
	#print('rank after symbol section', rank)
	
				
	#  2. Add repeat pattern sections (e.g. BCAB = [2,1,1], or BAAA = [3,1])
	#Note: FREEZE WARNING!  More symbols will combinatorically increase the iterations required, and can cause the algorithm to hang.
	sectionCount = countSymbolSection(k, n, s)
	for part in gen_int_parts_lifo(k,s):
		#print(part)
		if part == reps: break			
		rank += countSymbolReps(k, n, part)
	
	
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy, 
	#         so the rest is just one possible scheme to identify a sequence,
	#		  in an order that is similar to, but subtly different than lexical
	#  -------------
	
	# Calculate some values needed for ranking
	totalReps = len(reps)
	mspSeqTotal = multi_count(k, reps)
	mspRepsTotal = multi_count(totalReps, Counter(reps).values())
	combTotal = comb(totalSymbolsAvail, totalReps)
	repSectionSize = combTotal * mspSeqTotal * mspRepsTotal
	
	
	#  3. Add the combination section 
	combSectionSize = mspSeqTotal * mspRepsTotal
	combVals = getCombVals(totalSymbolsAvail, symbolsUsed)		
	combRank = vals_to_comb_rank(combVals) #COLEX order
	rank += (combSectionSize * combRank)	
	#print ('Comb Section size: ', combSectionSize)
	#print ('Comb Vals: ', combVals)
	#print('Comb Rank: ', combRank)
	
	
	#  4. Add the multi-set permutation rank of the repeats (mspReps)
	repPerm, repAlpha = getRepPerm(seq, totalSymbolsAvail, reps)	
	mspRepsRank = multi_perm_to_rank(repPerm, repAlpha)
	rank += (mspSeqTotal * mspRepsRank)	
	#print ('MSP Reps Perm: ', repPerm)
	#print ('MSP Reps Vals: ', repAlpha.vals)
	#print ('MSP Reps Rank: ', mspRepsRank)
	
	
	#  5. Add the multi-set permutation rank of the sequence (mspSeq)
	seqPerm, seqAlpha = seq_to_multi_perm(seq)
	mspSeqRank = multi_perm_to_rank(seqPerm, seqAlpha)
	rank += mspSeqRank
	#print ('MSP Seq Perm: ', seqPerm)
	#print ('MSP Seq Vals: ', seqAlpha.vals)
	#print ('MSP Seq Rank: ', mspSeqRank)	
	
	
	#print('Rank:', rank)
	return rank
	
	
def entropic_unrank(rank, totalSymbolsAvail, seqLen):	
	k = seqLen
	n = totalSymbolsAvail		
	
	#  1. Get symbol section
	#print('Starting Rank', rank)
	count = 0
	prevCount = count
	for i in range(1, totalSymbolsAvail+1):
		count += countSymbolSection(k, n, i)
		if count > rank:
			s = i			
			rank -= prevCount
			break
		else: prevCount = count
	try: s
	except: s = totalSymbolsAvail
	
	#print ('Symbol Section', s)
	#print ('Rank after symbol section:', rank)
				
	#  2. Get repeat pattern sections (e.g. BCAB = [2,1,1], or BAAA = [3,1])				
	sectionCount = 0
	prevCount = 0	
	for part in gen_int_parts_lifo(k, s):	
		#print('part', part)
		sectionCount += countSymbolReps(k, n, part)
		
		if sectionCount > rank:								
			reps = part	
			rank -= prevCount
			break			
		else: 			
			prevCount = sectionCount			
		
	#print ('Rank after rep section:', rank)
	
	#  -------------
	#   Note: Everything after this has the exact same amount of entropy	
	#  -------------
	
	# Now that we know the repeat pattern, we can calculate some values we'll need:
	totalReps = len(reps)
	mspRepsTotal = multi_count(totalReps, Counter(reps).values())
	mspSeqTotal = multi_count(k, reps)
	combTotal = comb(totalSymbolsAvail, totalReps)
	repSectionSize = combTotal * mspSeqTotal * mspRepsTotal
	#print('Reps', reps)
	#print ('mspSeqTotal', mspSeqTotal)
	#print ('mspRepsTotal', mspRepsTotal)
	#print ('combTotal', combTotal)
	#print ('Branch total:', mspSeqTotal * mspRepsTotal * combTotal)
	
	
	#  3. Get the combination from the combination rank 
	combSectionSize = repSectionSize // combTotal
	combRank = rank // combSectionSize 
	combVals = comb_rank_to_vals(combRank, totalSymbolsAvail, s) #COLEX order
	combVals.reverse()
	for i in range(len(combVals)):
		combVals[i] -= 1
	#print ('Comb Section size: ', combSectionSize)
	#print('Comb Rank: ', combRank)
	#print ('Comb Vals: ', combVals)
	
	
	#  4. Get the multi-set permutation of the repeats (mspReps)
	mspRepsSectionSize = combSectionSize // mspRepsTotal
	mspRepsRank = (rank % combSectionSize) // mspRepsSectionSize
	mspRepVals = getRepVals(reps)
	mspRepsPerm, mspRepsAlpha = multi_rank_to_perm(mspRepsRank, mspRepVals)	
	#print ('MSP Reps Vals: ', mspRepVals)
	#print ('MSP Reps Rank: ', mspRepsRank)
	#print ('MSP Reps Alpha: ', mspRepsAlpha.vals)
	#print ('MSP Reps Perm: ', mspRepsPerm)
	
	
	#  5. Add the multi-set permutation (MSP) rank of the repeat pattern	
	mspSeqRank = rank % mspRepsSectionSize 
	mspSeqVals = getSeqValsFromRepPerm(reps, mspRepsPerm, combVals)	
	mspSeqPerm, mspSeqAlpha = multi_rank_to_perm(mspSeqRank, mspSeqVals)
	seq = multi_perm_to_seq(mspSeqPerm, mspSeqAlpha)
	#print('MSP Seq Rank', mspSeqRank)
	#print('MSP Seq Vals', mspSeqVals)
	#print('MSP Seq Perm', mspSeqPerm)
	#print('MSP Seq Alpha', mspSeqAlpha.vals)
	#print (combRank, mspRepsRank, mspSeqRank)
	#print('Seq:', seq)
	return seq
	
#Combinatorial function to count the ways to rank:
# - Filling up a sequence of k length (e.g. k = 3)
# - Count of set of n symbols available (e.g. n = {A,B,C} = 3 )
# - Using subset s of symbols (e.g. s = 2)
def countSymbolSection(k, n, s): 
	if s < 0 or s > k or s > n: return 0	
	return factorial(s) * comb(n, s) * stirling2(k, s, exact=True)
	
	
#Combinatorial function to count ways to arrange repeat patterns with the given number of symbols
def countSymbolReps(k, totalSymbolsAvail, reps):
	totalReps = len(reps)
	
	mspSeq = multi_count(k, reps)
	mspReps = multi_count(totalReps, Counter(reps).values())
	combinations = comb(totalSymbolsAvail, totalReps)
	
	return combinations * mspSeq * mspReps
	

def getReps(v):
	s = set(v)	
	r = []
	for i in s:
		r.append(v.count(i))
	r.sort(reverse=True)
	return r

def getRepPerm(seq, symbols, reps):	
	repToCountId = {}	
	prev = reps[0]
	countId = 0
	
	for i in range(1, len(reps)):
		r = reps[i]
		
		if r != prev:
			repToCountId[prev] = countId					
			prev = r
			countId += 1
	
	repToCountId[prev] = countId		
	#print('rep to Count id', repToCountId)
	
	#Loop through symbols, and apply
	perm = []
	countsBySymbol = Counter(seq)
	for sym in range(symbols):
		if sym in countsBySymbol:
			rep = countsBySymbol[sym]
			countId = repToCountId[rep]			
			perm.append(countId)

	repVals = perm[:]
	repVals.sort() #vals are unordered, but sorted and stored in lexographic order for convenience		
	repAlpha = Alpha(repVals, 0)
	return perm, repAlpha	
	
def getRepVals(reps):		
	prev = reps[0]
	repVals = [0]
	repId = 0
	for i in range(1, len(reps)):
		r = reps[i]		
		if r != prev:
			repId += 1
			prev = r
			
		repVals.append(repId)
	return repVals	
	

def getSeqValsFromRepPerm(reps, mspRepsPerm, combVals):	
	#This function combines the symbols in the in the right quantities
	#to construct an alphabet that can then be used the the sequence
	#permutation to get the correct sequence back
	prev = reps[0]
	repIdToRepVal = [prev]	
	for i in range(1, len(reps)):
		r = reps[i]		
		if r != prev:
			repIdToRepVal.append(r)
			prev = r
							
	vals = []
	symId = 0	
	
	for repId in mspRepsPerm:
		rep = repIdToRepVal[repId]
		sym = combVals[symId]		
		vals += [sym] * rep
		symId += 1		
		
	return vals
	
def getCombVals(totalSymbolsAvail, symbolsUsed):
	vals = []
	for i in range(0, totalSymbolsAvail):		
		if i in symbolsUsed: vals.append(i+1)
	vals.reverse() #COLEX rank must be in descending order
	return vals