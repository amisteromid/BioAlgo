import numpy as np
#######################
spot={}
for i,a in enumerate('A C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y'.split()):
    spot[a]=i
matrix='''
4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
-2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
-1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
-2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
-2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
-1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
-1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
-1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
-1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
-2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
-1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
-1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
-1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
-3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
-2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7'''
matrix=list(map(lambda a:list(map(lambda x:int(x),a.split())),matrix.split('\n')[1:]))
########################


file = open("rosalind_ctea.txt").read()
file = file.split('>')
del file[0]
seqs=[]
for f in file:
    f = f.split('\n',1)[1].replace('\n','')
    seqs.append(f)
first_one = seqs[0]
second_one = seqs[1]


array=np.zeros( (len(first_one)+1 , len(second_one)+1) )
array[0,0]=np.NaN
array[0,1]=array[1,0]=-5
for i in range(2,len(second_one)+1):
    array[0,i]=array[0,i-1]
for i in range(2,len(first_one)+1):
    array[i,0]=array[i-1,0]

for a in range(1,len(first_one)+1):
    for d in range(1,len(second_one)+1):
        #print(a,d)
        if a==1 and d==1:
            array[a,d]=matrix[spot[first_one[a-1]]][spot[second_one[d-1]]]
        else:
            choose=[None for i in range(3)]
            choose[0]=array[a-1,d-1]+matrix[spot[first_one[a-1]]][spot[second_one[d-1]]]
            choose[1]=max(array[:a,d])-5
            choose[2]=max(array[a,:d])-5
            array[a,d]=max(choose)


print (array[len(first_one),len(second_one)])
