#Graph Coloring: 
import dwavebinarycsp
import networkx as nx
import matplotlib.pyplot as plt
import dimod
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

#list of ddepartments
dep = ['LP','OR','PT','CB','SC','BN','PA','TJ','CH']
#list of neighbors
neighbors = [('PA', 'LP'), ('PA', 'BN'), ('BN', 'LP'), ('BN','CB'), ('BN','SC'), ('CB','LP'), ('CB', 'OR'),('CB', 'PT'), ('CB', 'CH'),('CB','SC'),('OR','PT'),('PT', 'CH'),('PT', 'TJ'),('CH', 'TJ'),('SC','CH')]
#binary values for the colors
color_config = {(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0)}
colors =  len(color_config)


#graph graphic tools 
g = nx.Graph()
g.add_nodes_from(dep)
g.add_edges_from(neighbors)

map = {} #map to color

#used for making sure two neighboring departments
#don't have the same color 
def not_both(v, u):
    return not (v and u)

#plotting the map
def plot_map(sample):
	for d in dep: 
		for i in range(colors):
			if sample[d+str(i)]:
				map[d] = i #map to color
			
	node_colors = [map.get(node) for node in g.nodes()]
	nx.draw_circular(g, with_labels =True, node_color = node_colors, node_size=3000, cmap=plt.cm.rainbow)
	plt.show()

csp = dwavebinarycsp.ConstraintSatisfactionProblem(dwavebinarycsp.BINARY)

#single color 

for d in dep:
	var = [d+str(i) for i in range(colors)]
	csp.add_constraint(color_config,var)
# can't have the same color on neighboring departments
for neighbor in neighbors:
	v,u = neighbor
	for i in range(colors):
		var = [v+ str(i), u+str(i)]
		csp.add_constraint(not_both,var)
#add constraint : each pair of nodes with a shared edge can't have the same color
#for later....

bqm = dwavebinarycsp.stitch(csp)

sampler = EmbeddingComposite(DWaveSampler())         # doctest: +SKIP
response = sampler.sample(bqm, num_reads=80)         # doctest: +SKIP

sample = next(response.samples())      # doctest: +SKIP
if not csp.check(sample):              # doctest: +SKIP
    print("Failed to color map")
else:
    plot_map(sample)