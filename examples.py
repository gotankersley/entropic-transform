from entropic import entropic_rank, entropic_unrank
from math import log2

	
#https://sochima.me/set-shaping-theory
# compute the information content of a given string
def I(S):
	N = len(S)
	prob = [S.count(c)/N for c in S] # probability of each symbol in the string
	return round(-sum([log2(p) for p in prob]), 2)
	
TOTAL_SYM_AVAIL = 3
SEQ_LEN = 3
print('Small message examples:')
for r in range(0, TOTAL_SYM_AVAIL**SEQ_LEN):	
	seq = entropic_unrank(r, TOTAL_SYM_AVAIL, SEQ_LEN)
	rank = entropic_rank(seq, TOTAL_SYM_AVAIL)
	print (r, seq, rank, I(seq))	
	assert r == rank
	
	

TOTAL_SYM_AVAIL = 64
SEQ_LEN = 100
K_SHAPING = 3
print('Large example:')
print('Alphabet Symbols: ', TOTAL_SYM_AVAIL, '\n')

msg = [33, 51, 26, 44, 33, 61, 36, 28, 10, 48, 56, 45, 20, 36, 32, 43, 21, 28, 28, 23, 30, 3, 30, 51, 13, 37, 55, 34, 12, 36, 12, 47, 14, 25, 39, 31, 42, 46, 14, 47, 28, 4, 26, 26, 43, 43, 36, 47, 46, 52, 34, 57, 45, 2, 13, 35, 38, 45, 30, 42, 63, 29, 51, 57, 59, 20, 18, 33, 46, 53, 31, 28, 57, 19, 11, 60, 26, 63, 19, 4, 12, 59, 61, 16, 37, 3, 33, 2, 15, 19, 61, 35, 20, 52, 55, 33, 45, 25, 15, 15]
print('Random Message: ', msg, '\n')
print('Entropy of Random Message: ', I(msg), '\n')
#Convert msg into bits
msgBits = 2191682868323610206074412672256278751534977634292400501756475833529909753323059409647825439983699976059631680102414064582352842549313807925011638423879708831181008655333958010377167

entropicMsg = entropic_unrank(msgBits, TOTAL_SYM_AVAIL, SEQ_LEN + K_SHAPING)
print ('Unranked entropic msg:', entropicMsg)
print ('Entropy of msg: ', I(entropicMsg), '\n')
rank = entropic_rank(entropicMsg, TOTAL_SYM_AVAIL)
print ('Rank of sequence: ', rank)
assert msgBits == rank