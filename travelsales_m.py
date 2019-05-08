from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from itertools import combinations
Sampler = EmbeddingComposite(DWaveSampler())
#Eri Montano
#traveling sales man problem
#this is a 4 city problem and should work with more cities


#order of initial upper triangle for the city/weight order
distance = [4,6,5,1,4,7] 

E = len(distance) 

N = int((1+(1+E*8)**0.5)/2) #calcualte based on edges

cities = list(range(N))
 
gamma = 15.0 
chainstrength = 10
numruns = 100

edges = list(combinations(cities,2))

print (edges)


Q = {(a,b): 0 for a in range(E) for b in range(E)} 

for i in range(E): #form the diagonal
	Q[(i,i)] = distance[i]-(6*gamma)
	
for i in range(E):
	for j in range((i+1), E): #parse though the upper triangle
			row_edges = edges[i] 
			col_edges = edges[j]
			#the upper triangle will have some zeros and unless this condition meets..
			if row_edges[0] in col_edges or row_edges[1] in col_edges:
				#then we can calculate gamma 
				Q[(i,j)] = 2*gamma


print(Q) 


Response = EmbeddingComposite(DWaveSampler()).sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

## ------- Return results to user -------
'''
#R = iter(response)
#E = iter(response.data())
#for line in response:
#   sample = next(R)
#    S1 = [S[i] for i in sample if sample[i] > 0]
#    S0 = [S[i] for i in sample if sample[i] < 1]
#    print("S0 Sum: ", sum(S0), "\tS1 Sum: ", sum(S1), "\t", S0)
#Response = Sampler.sample_qubo(Q, num_reads=100)
'''
print(Response)
for index in range(len(Response.record.sample)):
    print("Sample:", Response.record.sample[index], "\tEnergy: ", Response.record.energy[index]) 