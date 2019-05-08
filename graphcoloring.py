#Graph Coloring: 
import dwavebinarycsp
import networkx as nx
import matplotlib.pyplot as plt
import dimod
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave_qbsolv import QBSolv


#gamma = 15.0 
chainstrength = 80
numruns = 20
cities = 5
#cities = {1,2,3,4,5,6,7,8,9}

na_cities = ['LP','OR','PT','CB','SC','BN','PA','TJ','CH']
na_neighbors = [('PA', 'LP'), ('PA', 'BN'), ('BN', 'LP'), ('BN','CB'), ('BN','SC'), ('CB','LP'), ('CB', 'OR'),('CB', 'PT'), ('CB', 'CH'),('CB','SC'),('OR','PT'),('PT', 'CH'),('PT', 'TJ'),('CH', 'TJ'),('SC','CH')]
color_config = {(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0)}

#neighbors = {'P':['E','A'], 'E':['A','D'], 'D':['S','A'], 'A':['S']}
neighbors = {1:[2,3,5], 2:[1,3,4], 3:[1,4,5], 5:[1,4]}
#neighbors = {1:[2,4,7,6], 4:[2,3,5,6], 2:[3], 5:[6,9], 3:[9,8], 8:[9]}
#neighbors = {1:[2,4,6,7], 2:[3,4], 3:[4,8,9], 4:[5,6], 5:[6,9], 6:[7], 8:[9]}
#--neighbors = {1:[2,4,7,6],2:[3,1,4], 3:[2,4,9,8], 4:[1,2,3,5,6], 5:[6,4,9], 6:[7,1,4,5], 7:[1,6], 8:[9,3], 9:[3,8,5]}


#neighbors = {0:[1,3,6,5],1:[2,0,3], 2:[1,3,8,7], 3:[0,1,2,4,5], 4:[5,3,8], 5:[6,0,3,4], 6:[0,5], 7:[8,2], 8:[2,7,4]}


g = nx.Graph()
g.add_nodes_from(na_cities)
g.add_edges_from(na_neighbors)

#colors =  len(color_config)

colors = 4
map = {}

N = cities * colors #number of ising variables

J = {} #coefficient of interaction
S = {} #adjacency matrix
h = [0] * N #magnetic field

A = 1

# Constraint Satisfaction Problem: requires all variables to have values
# This class adds all constraints and variables fefined fora problem and 
# provices funcitonality ti assust in the solution and ches if the solution 
# satisfies the constrainst (summary from dwave documentation) 

#csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)


for v in range(cities):
	for i in range(colors):
		h[v*colors+i] = -A
		
		for j in range(i+1, colors):
			if (v*colors+i, v*colors+j) not in J.keys():
				J[v*colors +i , v*colors +j] = 2*A
			else: 
				J[v*colors+i, v*colors+j] += 2*A
			
B = 2.5
for u, ne in neighbors.items():
	#print (neighbors.items())
	for i in range(colors):
		for v in (x-1 for x in ne):
			if u<v:
				for i in range(colors):
					if(u*colors+i, v*colors+i) not in J.keys():
						J[u*colors+i, v*colors+i] = B
					else:
						J[u*colors+i, v*colors+i] += B
						
# QUBO
Q = J.copy()
for i in range(N):
	if(i,i) not in Q.keys():
		Q[i,i] = h[i]
	Q[i,i] += h[i]
	
#adjacency matrix	
S= {}
for i in range(N):
	for j in range(N):
		if(i,j) in Q.keys() and Q[i,j] != 0:
			S[i,j]=1
			



response = QBSolv().sample_qubo(Q)
print("samples=" + str(list(response.samples())))
#print("energies=" + str(list(response.data_vectors['energy'])))


response1 = EmbeddingComposite(DWaveSampler()).sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)
print (response1) 


##thing
for index in range(len(response1.record.sample)):
    print("Sample:", response1.record.sample[index], "\tEnergy: ", response1.record.energy[index]) 
#end of thing


sample = next(response1.samples())

# plot map

#for c in range(cities): 
#	for i in range(colors):
#		if sample[c + i]:
#			map[c] = i 
	
#	node_colors = [map.get(node) for node in g.nodes()]
#	nx.draw_circular(g, with_labels =True, node_color = node_colors, node_size=3000, cmap=plt.cm.rainbow)
#	plt.show()


# format respomse dwave : matrix 
