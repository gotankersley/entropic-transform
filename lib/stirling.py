from math import factorial, comb
from functools import lru_cache
from functools import cmp_to_key #https://wiki.python.org/moin/HowTo/Sorting/

#Set Part canonical order: First element of each set, (each set is already assumed to be ordered within itself)
def setPartCompare(a, b):	
	if a[0] > b[0]: return 1
	elif a[0] < b[0]: return -1
	else: return 0	
	
	
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
	

def seq_to_set_part(symSeq, symCount):
	setPart = []	
	for i in range(symCount):
		setPart.append([])
			
	for pos in range(len(symSeq)):
		sym = symSeq[pos]
		setPart[sym].append(pos)
		
	#Sort in canonical order
	return sorted(setPart, key=cmp_to_key(setPartCompare))		

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

#stirCache = {}
@lru_cache(None)
def stirling2MaxLessThan(n, k, m):
	#global stirCache
	#stirKey = str(n) + ',' + str(k) + ','+ str(m)
	#if stirKey in stirCache: return stirCache[stirKey]
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
	#stirCache[stirKey] = S[n][k]
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
	
def getSizeOfLargestPart(setPart):
	largestPart = -1
	for p in setPart:
		largestPart = max(largestPart, len(p))
	return largestPart