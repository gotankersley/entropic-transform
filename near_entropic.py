from scipy.special import stirling2
from math import comb, factorial

from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/

from lib.myrvold import *
from lib.combinations import *
from lib.stirling import *



#ABOUT:  This is a function to rank a sequence in NEAR sorted entropic order, rather than the regular place-value base system.
# This means that sequences with the less entropy, (i.e. fewer distinct symbols), like 'AAA', 'BBB', are ranked first,
# and those with the more entropy, (i.e. more differing symbols), like 'ABCD', are ranked last.
#
# Sequences are uniquely identified in the following way:
#  1. By symbol section, (i.e. all the sequences with only 1 symbol are grouped together, and so on)
#  2. By Set Partition rank (Stirling number of the second kind), in the lexicographical order 
#	  given by the Restrictive Growth Function (RGF).   
#  -------------
#   Note 1: The lexicographical order for set Partitions is somewhat different than the 
#		  integer partition order needed for ranking in exact entropic order, which is why
#		  this can only be considered NEAR sorted entropic order. 
#		  For example:	CCDD (entropy 4) comes BEFORE BABB (entropy 3.25)
#		  However, it IS sorted by section, so everything in a preceeding section has
#		  less entropy than all the sections following, so it is correctly ordered that way.
#	
#	Note 2: Everything after this has the exact same amount of entropy, 
#         so the rest is just one possible scheme to identify the sequence
#  -------------
#  3. By the combination rank of the symbols used from the set of total available symbols
#  4. By the permutation rank of the mapping the symbols to the set partition
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
	

#Combinatorial function to count the ways to rank:
# - Filling up a sequence of k length (e.g. k = 3)
# - Count of set of n symbols available (e.g. n = {A,B,C} = 3 )
# - Using subset s of symbols (e.g. s = 2)
def countSymbolSection(k, n, s): 
	if s < 0 or s > k or s > n: return 0	
	return comb(n, s) * factorial(s) * stirling2(k, s, exact=True)
	
	

def getCombVals(vals):
	#COLEX rank must be in descending order
	combVals = []
	for i in range(len(vals)-1, -1, -1):
		combVals.append(vals[i]+1)
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