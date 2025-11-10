from near_entropic import near_entropic_rank, near_entropic_unrank
from lib.baselib import b2n
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
	seq = near_entropic_unrank(r, TOTAL_SYM_AVAIL, SEQ_LEN)
	rank = near_entropic_rank(seq, TOTAL_SYM_AVAIL)
	print (r, seq, rank, I(seq))	
	assert r == rank
	
print ('---')	

msg = [33, 51, 26, 44, 33, 61, 36, 28, 10, 48, 56, 45, 20, 36, 32, 43, 21, 28, 28, 23, 30, 3, 30, 51, 13, 37, 55, 34, 12, 36, 12, 47, 14, 25, 39, 31, 42, 46, 14, 47, 28, 4, 26, 26, 43, 43, 36, 47, 46, 52, 34, 57, 45, 2, 13, 35, 38, 45, 30, 42, 63, 29, 51, 57, 59, 20, 18, 33, 46, 53, 31, 28, 57, 19, 11, 60, 26, 63, 19, 4, 12, 59, 61, 16, 37, 3, 33, 2, 15, 19, 61, 35, 20, 52, 55, 33, 45, 25, 15, 15]
TOTAL_SYM_AVAIL = 64
SEQ_LEN = len(msg)
K_SHAPING = 1

print('Random Message: ', msg)
print('Alphabet Symbols: ', TOTAL_SYM_AVAIL)
print('Entropy of Random Message: ', I(msg), '\r\n')
msgRank = b2n(TOTAL_SYM_AVAIL, msg)

print ('---')	

print('Starting from Base-A rank:')
entropicMsg = near_entropic_unrank(msgRank, TOTAL_SYM_AVAIL, SEQ_LEN + K_SHAPING)
print ('Unranked entropic msg:', entropicMsg)
print ('Entropy of msg: ', I(entropicMsg), '\n')
rank = near_entropic_rank(entropicMsg, TOTAL_SYM_AVAIL)

assert msgRank == rank

print ('---')

print('Starting from Entropic rank:')
rank2 = near_entropic_rank(msg, TOTAL_SYM_AVAIL)
shapedSeq = near_entropic_unrank(rank2, TOTAL_SYM_AVAIL, SEQ_LEN + K_SHAPING)
print ('Unranked shaped seq:', shapedSeq)
print ('Entropy of shaped seq: ', I(shapedSeq))