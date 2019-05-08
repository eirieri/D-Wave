# Exercise 1:  Box problem
# Find the minimum sum of only two numbers from the list
## ------- Set up our list of numbers -------
S = [4, 21, 3, 17]

## ------- Run our QUBO on the QPU -------
# Set up QPU parameters
chainstrength = 50
numruns = 10
gam = 30

Q = {}
#diagonal of the matrix 
for i in range(len(S)):
	Q[(i,i)] = (S[i]-3*gam)
	
#upper triangle of the matrix
for i in range(len(S)):
	for j in range((i+1), len(S)):
			Q[(i,j)] = 2*gam


# Run the QUBO on the solver from your config file
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
response = EmbeddingComposite(DWaveSampler()).sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

## ------- Return results to user -------
R = iter(response)
E = iter(response.data())
for line in response:
    sample = next(R)
    S1 = [S[i] for i in sample if sample[i] > 0]
    S0 = [S[i] for i in sample if sample[i] < 1]
    print("S1 Sum: ", sum(S1), "\t", S1)
