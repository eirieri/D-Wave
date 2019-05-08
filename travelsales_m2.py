from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from itertools import combinations
Sampler = EmbeddingComposite(DWaveSampler())


#traveling salesman problem 

distance = [4,6,5,1,4,7]

E = len(distance)

N = int((1+(1+E*8)**0.5)/2) #calcualte based on edges

cities = list(range(N))
 
gamma = 15.0 
chainstrength = 10
numruns = 100

edges = list(combinations(cities,2))

print (edges)

#x_ij=1 if city i is in the jth stop
#path to minimize, for each pair (from stop 1 to stop2) : sum_ij(x_i,1*xj,2*d_ij)... sum for all pairs->d_ij
#constraint: each city is visited only once: 
#gamma*(sum_ij x_ij -1)^2= gamma*sum_ij (x_ij^2)-gamma*2*sum_ij(x_ij)+sum_ij 1-> gamma -2gamma= -gamma
#constraint: only one visit in each stop: 
#->-gamma


Q = {(a,b): 0 for a in range(E) for b in range(E)} #init every pair to zero

for i in range(E):
	Q[(i,i)] = distance[i]-(6*gamma)
	
for i in range(E):
	for j in range((i+1), E):
			row_edges = edges[i]
			col_edges = edges[j]
			if row_edges[0] in col_edges or row_edges[1] in col_edges:
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