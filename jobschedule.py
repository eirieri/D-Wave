from sympy import *
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

init_printing(use_latex=True)

# Functions to make QUBO matrix using sympy
def monomToIndicies(monom):
    # Only one power
    if sum(monom) == 1:
        for i in range(len(monom)):
            if monom[i] == 1:
                return i, i
    
    
    # There must be a quadratic monom
    elif sum(monom) == 2:
        # Find the first non-zero index
        for i in range(len(monom)):
            if monom[i] == 2:
                return i, i
            elif monom[i] == 1:
                break

        # Find the second non-zero index
        for j in range(i + 1, len(monom)):
            if monom[j] == 1:
                return i, j
    
    else:
        return None, None


def makeCoeffMatrix(poly, n):
    coeffs = poly.coeffs()
    monoms = poly.monoms()

    if len(coeffs) != len(monoms):
        print("Can't make matrix")
        return None

#     print([sum(m) for m in monoms])

    # Init nxn matrix
    m = [
        [0]*n for _ in range(n)
    ]

    for i in range(len(monoms)):
        x, y = monomToIndicies(monoms[i])
        if x is None:
            continue
        m[x][y] += coeffs[i]

    return Matrix(m)


# Function to parse output of DWave solver
def parse(sample, m, N, W, L):
    for a in range(m):
        s = f"\tMachine {a}: "
        total = 0
        for i in range(N):
            if sample[a*N + i] == 1:
                s += f"{i}, "
                total += L[i]
        s += f"\tTotal = {total}"
        print(s)
    
    for a in range(1, m):
        total = 0
        for i in range(W + 1):
            if sample[m*N + (a - 1)*(W + 1) + i] == 1:
                total += 2**i
        
        s = f"\tY Matrix: {total}"
        print(s)


##############################################################################
# Change Stuff Here
##############################################################################
# Looping symbols for sigmas
i, n, a = symbols('i n a', cls=Idx)
# Vector symbols
L = symbols('L', cls=IndexedBase)
z = symbols('z', cls=IndexedBase)  # For x matrix
w = symbols('w', cls=IndexedBase)  # For y matrix

chainstrength = 750
numruns = 100
gamma1 = 70
gamma2 = 10
B = 1

# List of job lengths
# LL = [2, 9, 3, 3, 6, 7]
LL = [3, 4, 6]


N = len(LL)  # Number of jobs
m = 2  # Number of machines
fancy_m = N * max(LL)  # Worse case for y matrix
# This value should be the highest power of 2 such that 2**W < fancy_m
W = 4


# Constraint for every job must be on one and only one machine
const1 = expand(summation((2*summation(z[N*a + i], (a, 0, m - 1)) - 1)**2, (i, 0, N - 1)))
# Constraint to make machine 0 the longest
const2 = expand(summation((summation(2**n*w[(W + 1)*a + n], (n, 0, W - 1)) + (fancy_m + 1 - 2**W)*w[(W + 1)*a + W] + summation(L[i]*(z[N*a + i] - z[i]), (i, 0, N - 1)))**2, (a, 1, m - 1)))

# Objective function: Minimize machine 0's length
obj1 = expand(summation(L[i]*z[i], (i, 0, N - 1)))

# Combine contraints and objective functions into one polynomial
poly = Poly(gamma1*const1 + gamma2*const2 + B*obj1, [z[i] for i in range(N*m)] + [w[i] for i in range(W + 1, m*(W + 1))])
mat = makeCoeffMatrix(poly, m*N + (m - 1)*(W + 1))

subs = {}
for ii in range(len(LL)):
    subs[L[ii]] = LL[ii]

mat_subs = mat.subs(subs)

# Move matrix into Q dictionary for DWave solver.
Q = {}
for ii in range(m*N + (m - 1)*(W + 1)):
    for j in range(m*N + (m - 1)*(W + 1)):
        Q[(ii, j)] = mat_subs[ii, j]

# Solve
response = EmbeddingComposite(DWaveSampler()).sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

# Parse output of DWave
print(response)
for index in range(len(response.record.sample)):
    sample = response.record.sample[index]
    print("Sample:", sample, "\tEnergy:", response.record.energy[index])
    if (response.record.energy[index]) < 100:
        parse(sample, m, N, W, LL)
