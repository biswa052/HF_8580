import numpy as np
import math

def kinetic(basis, S):
	T = np.zeros((4,4))
	i = 0
	j = 0
	q = 0.
	Q = 0.
	p = 0.
	while i < 4:
		while j < 4:
			Q = basis[i].z - basis[j].z
			Nind1 = 0
			while Nind1 < len(basis[i].alpha):
				Nind2 = 0
				while Nind2 < len(basis[j].alpha):
					p = basis[i].alpha[Nind1] + basis[j].alpha[Nind2]
					q = basis[i].alpha[Nind1] * basis[j].alpha[Nind2]/p
					T[i,j] += q*(3.-(2.*q*Q*Q))*(basis[i].N[Nind1]*basis[j].N[Nind2]*((math.pi/p)**(1.5))*(math.exp(-q*Q*Q))*basis[i].c[Nind1]*basis[j].c[Nind2])
					Nind2 += 1
				Nind1 += 1
			j += 1
		j = 0
		i += 1
	
	return T