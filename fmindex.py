import numpy as np
import random

def getUnique(L):
    newL = []
    for l in L:
        if(l not in newL): newL.append(l)
    return newL

def bwt_permute(S):
    s = S
    s_list = []
    s_list.append(s)
    
    for _ in range(0, len(s)-1):
        s = s[1:] + s[0]
        s_list.append(s)    
    sa_t = sorted(range(len(s_list)), key=lambda k: s_list[k])
    s_list.sort()    
    L = ''
    F = ''
    for i in range(0, len(s_list)):
        L += s_list[i][-1]
        F += s_list[i][0]
    return s_list, sa_t, F, L, s_list.index(S)

def bwt_unpermute(S, l_idx):    
    s_list = [ x for x in S ]
    s_list.sort()
    for _ in range(0, len(S)-1):
        for j in range(0, len(S)):
            s_list[j] = S[j] + s_list[j]
        s_list.sort()
    return s_list[l_idx]

def getC(x):
    dic = {}
    ex_c = ''
    for cnt, c in enumerate(x):
        if(c != ex_c):
            dic[c] = cnt
            ex_c = c
    return dic

def occ(L):
    unique_c = list(set(L))
    unique_c.sort()
    dic = dict()    
    for c in unique_c:
        dic[c] = np.array([0]*len(L))
    for idx,l in enumerate(L):        
        dic[l][idx] += 1
        if(idx < len(L)-1):
            for c in unique_c:
                dic[c][idx+1]  = dic[c][idx]    
    return dic

def generateT(length):
    t = ''
    for _ in range(0,length-1):
        t += random.sample(['A','G','C','T'],1)[0]
    return t+'$'

def generateP(t,length):
    sp = random.randint(0,len(t)-length)
    return t[sp:(sp+length)]

def printm(m):
    for idx,x in enumerate(m):
        print(idx,':',x)

def calculateD(W,X,C):
    m,_,f,B2,_ = bwt_permute(X[::-1])
    #printm(m)
    O2 = occ(B2)
    #print(O2)
    
    D = []
    k = 1
    l = len(X)-1
    z = 0
    
    flg = True
    for w in W:
        if(w not in C):
            z+=1
            D.append(z)
            flg = True
            l = len(W)-1
            continue
            
        if(flg):
            k = C[w] + 1 - 1
            flg = False
        else:
            k = C[w] + O2[w][k-1] + 1 - 1
        l = C[w] + O2[w][l]        
        #print(w, k,l)        
        if(k > l):
            z+=1
            D.append(z)
            flg = True
            l = len(W)-1
        else:
            D.append(z)
    
    return D

def inexRecur(W,i,z,k,l,C,O,D,flg):
        
    if( z<0 or z<D[i] ): return None    
    if( i<0 ): 
        print('Catch It!')
        return [[k,l]]
    
    I = []
    
    for b in ['A','C','G','T']:
        if(flg):
            k = C[b]-1 + 1
        else:
            k = C[b]-1 + O[b][k-1] + 1
        l = C[b] + O[b][l]
        
        #print('b',b,'k',k,'l',l)
        if(k <= l):
            if(b==W[i]):
                print(i,b,W[i],'Match','k:',k,'l:',l)
                print(W,i-1,z,k,l,D,False)
                ret = inexRecur(W,i-1,z,k,l,C,O,D,False)
                if(ret is not None):
                    print(ret)
                    I += ret
            else:
                #print(i,b,W[i],'Different')
                ret = inexRecur(W,i-1,z-1,k,l,C,O,D,False)
                if(ret is not None):
                    print(ret)                    
                    I += ret
                    
    return (None if len(I)==0 else I)

def inexactSearch(X,W,z,B,C,O):
    D = calculateD(W,X,C)    
    #print(D)
    return inexRecur(W,len(W)-1,z,1,len(X)-1, C, O, D,True)


#t = generateT(10)
#p = generateP(t,4)

#BWT from reference.
X = generateT(11)
'TTGCGACAGG$'
m, sa_t, f, B, l_idx = bwt_permute(X)
c = getC(f)
o = occ(B)

W = 'GATA'

result = inexactSearch(X,W,1,B,c,o)
result = getUnique(result)
printm(m)
