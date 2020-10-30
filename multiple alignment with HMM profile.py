import numpy as np
from numpy import argmax,log
import pandas as pd

####################################
##construct a profile HMM from sequences##
####################################
def hmmbuilder(theta,sigma,x,seqs):
    seqs=list(map(lambda a:[i for i in a],seqs))
    align = list(zip(*(seqs)))
    Ms=[]
    for i,a in enumerate(align):
        if a.count('-')/len(seqs)<theta:
            Ms.append(i)

    
    mstate =['M{} D{} I{}'.format(i,i,i).split() for i in range(1,len(Ms)+1)]
    states = ['S', 'I0'] + [each for m in mstate for each in m] + ['E']
    
    transitions = pd.DataFrame(data = 0.0 , index=states , columns=states)
    emissions = pd.DataFrame(data = 0.0 , index=states , columns=x)

    sray=np.array(seqs,dtype='<U4')
    counter=0
    for i in range(len(sray[0])):
        if i in Ms:
            counter += 1
            for ii,seq in enumerate(seqs):
                if seq[i]=='-' :
                    sray[ii,i]= 'D'+str(counter)
                else:
                    sray[ii,i]='M'+str(counter)
        else:
            for ii,seq in enumerate(seqs):
                if seq[i] != '-': sray[ii,i]='I'+str(counter)
            
    for ii,seq in enumerate(sray):
        skip=0
        for i in range(len(seq)+1):
            if i<len(seq) and seq[i]=='-':
                skip +=1
            else: 
                if i==0 or skip==i:
                    transitions.loc['S',seq[i]] += 1
                    if seq[i][0] != 'D': emissions.loc[seq[i],seqs[ii][i]] += 1
                elif i==len(seq):
                    transitions.loc[seq[i-1-skip],'E'] += 1
                else:
                    if seq[i][0] != 'D': emissions.loc[seq[i],seqs[ii][i]] += 1
                    transitions.loc[seq[i-1-skip],seq[i]] += 1
                skip=0
    
    transitions = transitions.div(transitions.sum(1) + 1e-10, axis=0).round(3)
    emissions = emissions.div(emissions.sum(1) + 1e-10, axis=0).round(3)

    fillit={'S':['I0', 'M1','D1'],'I0':['I0','M1','D1'],'E':[]}
    for s in states:
        if s not in fillit:
            if s[1:]==str(len(Ms)):
                fillit[s]=['I'+str(s[1:]),'E']
            else:
                if s[0]=='I':
                    fillit[s]=[s,'M'+  str(int(s[1:])+1) ,'D'+str(int(s[1:])+1)]
                else:
                    fillit[s]=['I'+str(s[1:]),'M'+str(int(s[1:])+1),'D'+str(int(s[1:])+1)]
    for s in states:
        for ss in fillit[s]:
            transitions.loc[s,ss] += sigma
    for s in states:
        if s[0]=='M' or s[0]=='I':
            for o in x:
                emissions.loc[s,o] += sigma
    transitions = transitions.div(transitions.sum(1) + 1e-10, axis=0) + 1e-100
    emissions = emissions.div(emissions.sum(1) + 1e-10, axis=0)  + 1e-100
    
    #transitions.to_csv('trans matrix.csv')
    #emissions.to_csv('emis matrix.csv')
    
    return transitions,emissions,states,x

####################################
###Align a sequence with a HMM profile####
####################################
def hmmaligner(transitions,emissions,pi,x,wanted):
    matches = int((len(transitions)-3)/3)
    
    dic={(state,i):[-np.inf,''] for state in pi for i in range(len(wanted)+1)}
    seq = [i for i in range(len(wanted)+1)]
    scoring = pd.DataFrame( data = 'float_tuple' , index = seq , columns = pi)
    #scoring.loc[0,'S']=(-np.inf,['M',1])

    for i in scoring.columns[:]:
        if i[0] != 'D': scoring.loc[0,i]='XXXXXXXXXXXXXX'

    
    # from S #####################
    scoring.loc[1,'M1']=( log(transitions.loc['S','M1'] )+ log(emissions.loc['M1',wanted[0]]),(0,'S'))
    scoring.loc[1,'I0']=( log(transitions.loc['S','I0'] )+ log(emissions.loc['I0',wanted[0]]),(0,'S'))
    scoring.loc[0,'D1']=( log(transitions.loc['S','D1'] ),(0,'S'))
    #initialize I0,M1,D1 row ################
    for i in range(2,len(wanted)+1):
        scoring.loc[i,'I0']=( scoring.loc[i-1,'I0'][0]+log(transitions.loc['I0','I0'] )+ log(emissions.loc['I0',wanted[i-1]]),(i-1,'I0'))
        scoring.loc[i,'M1']=( scoring.loc[i-1,'I0'][0]+log(transitions.loc['I0','M1'] )+ log(emissions.loc['M1',wanted[i-1]]),(i-1,'I0'))
    for i in range(2,len(wanted)+2):
        scoring.loc[i-1,'D1']=( scoring.loc[i-1,'I0'][0]+log(transitions.loc['I0','D1'] ),(i-1,'I0'))
    #initialize 0 column just D is allowed #############
    for i in range(2,matches+1):
        scoring.loc[0,'D'+str(i)]= (scoring.loc[0,'D'+str(i-1)][0]+log(transitions.loc['D'+str(i-1),'D'+str(i)]),(0,'D'+str(i-1))) 
    #initialize 1 column of I ####################
    for i in range(1,matches+1):
        scoring.loc[1,'I'+str(i)] = (scoring.loc[0,'D'+str(i)][0] + log(transitions.loc['D'+str(i),'I'+str(i)])+log(emissions.loc['I'+str(i),wanted[i-1]]),(0,'D'+str(i)))
    #initialize 1 column of M ####################
    for i in range(2,matches+1):
        scoring.loc[1,'M'+str(i)] = (scoring.loc[0,'D'+str(i-1)][0] + log(transitions.loc['D'+str(i-1),'M'+str(i)])+log(emissions.loc['M'+str(i),wanted[i-1]]),(0,'D'+str(i-1)))

    def recurrence_M(l,k):
        assert l>1 and k>1
        d_prev_id = 'D'+str(l-1)
        m_prev_id = 'M'+str(l-1)
        i_prev_id = 'I'+str(l-1)
        m_curr_id = 'M'+str(l)
        k_id = wanted[k-1]
        score_M = scoring.loc[k-1,'M'+str(l-1)][0] + log(emissions.loc[m_curr_id,k_id]*transitions.loc[m_prev_id,m_curr_id])
        score_D = scoring.loc[k-1,'D'+str(l-1)][0] + log(emissions.loc[m_curr_id,k_id]*transitions.loc[d_prev_id,m_curr_id])
        score_I = scoring.loc[k-1,'I'+str(l-1)][0] + log(emissions.loc[m_curr_id,k_id]*transitions.loc[i_prev_id,m_curr_id])
        score_max = max(score_M,score_D,score_I)
        if score_M == score_max:
            return (score_M,(k-1,'M'+str(l-1)))
        if score_D == score_max:
            return (score_D,(k-1,'D'+str(l-1)))
        if score_I == score_max:
            return (score_I,(k-1,'I'+str(l-1)))
    
    def recurrence_D(l,k):
        d_prev_id = 'D'+str(l-1)
        m_prev_id = 'M'+str(l-1)
        i_prev_id = 'I'+str(l-1)
        d_curr_id = 'D'+str(l)
        score_M = scoring.loc[k,'M'+str(l-1)][0] + log(transitions.loc[m_prev_id,d_curr_id])
        score_D = scoring.loc[k,'D'+str(l-1)][0] +log(transitions.loc[d_prev_id,d_curr_id])
        score_I = scoring.loc[k,'I'+str(l-1)][0] + log(transitions.loc[i_prev_id,d_curr_id])
        score_max = max(score_M,score_D,score_I)
        if score_M == score_max:
            return (score_M,(k,'M'+str(l-1)))
        if score_D == score_max:
            return (score_D,(k,'D'+str(l-1)))
        if score_I == score_max:
            return (score_I,(k,'I'+str(l-1)))

    def recurrence_I(l,k):
        d_curr_id = 'D'+str(l)
        m_curr_id = 'M'+str(l)
        i_curr_id = 'I'+str(l)
        k_id = wanted[k-1]
        score_M = scoring.loc[k-1,'M'+str(l)][0] + log(emissions.loc[i_curr_id,k_id]*transitions.loc[m_curr_id,i_curr_id])
        score_D = scoring.loc[k-1,'D'+str(l)][0] + log(emissions.loc[i_curr_id,k_id]*transitions.loc[d_curr_id,i_curr_id])
        score_I = scoring.loc[k-1,'I'+str(l)][0] +\
                  log(emissions.loc[i_curr_id,k_id]*transitions.loc[i_curr_id,i_curr_id])
        score_max = max(score_M,score_D,score_I)
        if score_M == score_max:
            return (score_M,(k-1,'M'+str(l)))
        if score_D == score_max:
            return (score_D,(k-1,'D'+str(l)))
        if score_I == score_max:
            return (score_I,(k-1,'I'+str(l))) 
    #___________fill it up_____________________________

    scoring.to_csv('scoring.csv')    
    for i in range(2,len(wanted)+1):
        scoring.loc[i,'I1']=recurrence_I(1,i)
    
    for i in range(2,matches+1):
        scoring.loc[1,'D'+str(i)]= recurrence_D(i,1)

    for l in range(2,matches+1):
        for k in range(2,len(wanted)+1):
            scoring.loc[k,'M'+str(l)] = recurrence_M(l,k)
            scoring.loc[k,'D'+str(l)] = recurrence_D(l,k)
            scoring.loc[k,'I'+str(l)] = recurrence_I(l,k)

    #scoring.to_csv('scoring.csv')
            
    #_________________backtrack_____________________
    final_score=[]
    for x in ['M','D','I']:
        final_score.append(scoring.loc[len(wanted),x+str(matches)][0]+log(transitions.loc[x+str(matches),'E']))
    last_state=(['M','D','I'][argmax(final_score)]+str(matches))
    last_state=(len(wanted),last_state)

    path=[]
    while last_state[1] != 'S':
        path.append(last_state[1])
        last_state = scoring.loc[last_state][1]
    print (' '.join([i for i in reversed(path)]))
        



###############
f = open('rosalind_ba10g.txt').read().split('\n',6)
wanted = f[0].rstrip()
theta, sigma = [float(i) for i in f[2].split()]
x = f[4].rstrip().split()
seqs = f[6].rstrip().split()
hmmaligner(*(hmmbuilder(theta,sigma,x,seqs)),wanted)
    
