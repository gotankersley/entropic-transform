from collections import deque 

def gen_int_parts_lifo(m, n):	
	parts = [m-(n-1)] + [1]*(n-1)	
	stack = deque()
	
	stack.append(parts)
	used = set()
	while stack:
		
		top = stack.popleft()		
		yield top[:]
		for i in range (0, n):
			if top[i] <= 2: break
			for j in range(i+1, n):				
				if top[i]-1 >= top[j]+1:
					if top[j]+1 <= top[j-1]:
						if top[i]-1 >= top[i+1]:
							topCopy = top[:]
							topCopy[i] -= 1
							topCopy[j] += 1
							topCopyStr = str(topCopy)
							if topCopyStr not in used:								
								stack.append(topCopy)
								used.add(topCopyStr)
					if top[j] == 1: break