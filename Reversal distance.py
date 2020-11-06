file=open("rosalind_rear.txt").read()
file=file.split('\n\n')
pairs=[]

for f in file:
    f=f.split('\n')
    if f.count('')==1: del f[-1]
    pair=[]
    for each in f:
        each=each.split(' ')
        pair.append(tuple(map(int,each)) )
    pairs.append(pair)
    
def best_intersection(P, Q):
  for perm in P:
    if perm in Q:
      return P[perm] + Q[perm]
  return -1

def make_step(processed, border, next_border):
  #print("processed : ",processed)
  #print("border ; ",border)
  #print("next_border : ",next_border)
  for perm in border:
    for j in range(len(perm)+1):
      for i in range(j-1): # at least 2 numbers reversed
        new_perm = perm[:i] + perm[i:j][::-1] + perm[j:]  # reverse [i,j]

        if new_perm not in processed and new_perm not in border:
          next_border[new_perm] = border[perm] + 1

  processed.update(border)
  border.clear()
  border.update(next_border)
  next_border.clear()
  #print("processed : ",processed)
  #print("border ; ",border)
  #print("next_border : ",next_border)



def rev_dist(perm1, perm2):
  processed1 = {}
  border1 = {}
  next_border1 = {}

  processed2 = {}
  border2 = {}
  next_border2 = {}

  border1[tuple(perm1)] = 0
  border2[tuple(perm2)] = 0

  lamp = False
  while best_intersection(border1, border2) == -1:
    if lamp:
      #print ("lampO")
      make_step(processed1, border1, next_border1)
    else:
      #print ("lampf")
      make_step(processed2, border2, next_border2)
    lamp = not lamp

  return best_intersection(border1, border2)



for pair in pairs:
    print(rev_dist(pair[0],pair[1]),end=' ')
