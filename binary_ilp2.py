#Bianry Integer Linear Programming:

#Let x1, x2 ... xn be N binary variables turn to vector x
# vector x is our solution as binary numbers 0-1 

import numpy as np
from sympy import *
from dwave_qbsolv import QBSolv

## ------- Run QUBO on the QPU -------
# Set up QPU parameters
chainstrength = 80
numruns = 100
gam = 200
gam2 = 300
gam3 = 80

#HARD CODED PROBLEM:
#initializing symbols, array of 4 for vector x
x1, x2, x3 ,x4 = symbols('x1 x2 x3 x4')
#given matrix s: 
#s = Matrix([[58,44,26,23],[25,29,13,17],[43,25,23,29]])
s =np.matrix([[58,44,26,23],[25,29,13,17],[43,25,23,29]])

#expected return value also given as array c

c = [217,125,88,109]
x = Matrix([x1,x2,x3,x4])

#goal: maximize the return value, c is multiplied by x
#	   to obtain the result for vector x
print ('vector max: c*x')
max = np.dot(c,x)
print (max)

#constraint for limited funds of each month
f = [107,55,95]

#creating matrix b when b=s*x
print ('vector b: s*x')
b = np.dot(s,x)
print (b)

#i=0
#b1 = s.sum(axis=1)
#print(b1)

#all zero bc project 2 need to be kecked out

#for i in range(len(b1)):
#	if b1[i] < f[i]:
#		print(f[i])
#		x[i] = 1
#	else:
#		x[i]=0
#	print (b1[i])
#	print (x[i])

	
#still need to parse it in rows

#QUBO matrix

Q = {}
for i in range(len(s)):
	Q[(i,i)] = (c[i]-587*gam2)
	#Q[(i,i)] = (c[i]-1075*gam2) + c[i]
#Adjacency Matrix
for i in range(len(s)):
	for j in range((i+1), len(c)):
			Q[(i,j)] = (2*gam3)


#attempt:2
#Q = {}
#for row in s:

#	for i in range(len(c)):
#		Q[(i,i)] = s[i]-587*gam
#		print (c[i])
	
#	for i in range(len(c)):	
#		for j in range((i+1), len(c)):
#			Q[(i,j)] = 6*gam


from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
response = EmbeddingComposite(DWaveSampler()).sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

## ------- Return results to user -------
R = iter(response)
print(response)



resp = QBSolv().sample_qubo(Q)
print("samples=" + str(list(resp.samples())))
#print("energies=" + str(list(resp.data_vectors['energy'])))


E = iter(response.data())
#print(E)

for line in response:
#second loop to go through each row 
    sample = next(R)
    S1 = [c[i] for i in sample if sample[i] > 0]
    S0 = [c[i] for i in sample if sample[i] < 1]
    print("S1 Sum: ", sum(S1), "\t", S1)
	#need to print in binary 
	
