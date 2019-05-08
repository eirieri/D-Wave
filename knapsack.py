## Knapsack with integer weights
## ------- Set up our list of numbers -------
import math

#list of objects
S = [10,3,5,7]
#value of objects
C = [6,7,2,4]
#lenght of list of object 
N = len(S)
#capacity of knapsack
W = 15

#should need for QUBO but it's used in part 2
spins = N + (1+ int(math.log10(W)))

#CONSTRAINTS

## ------- Run our QUBO on the QPU -------
# Set up QPU parameters
chainstrength = 1500
numruns = 40
gam = 100

Q = {}
for i in range(len(S)):
	Q[(i,i)] = (S[i]*S[i]-30*S[i])*gam - C[i]
	
for i in range(len(S)):
	for j in range((i+1), len(S)):
			Q[(i,j)] = (2*S[i]*S[j]+2)*gam


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
