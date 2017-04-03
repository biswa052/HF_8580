import numpy as np
import math 

def F0(t):
	if t == 0:
		F0ans = 1
	else:
		F0ans = .5 * (math.pi/t)**(.5) * math.erf(t**.5)
	return F0ans

def attract(basis, coords):
	V_eN = np.zeros((4,4))
	i = 0
	j = 0
	q = 0.
	Q = 0.
	p = 0.
	P = 0.
	zH = coords.item((1,2))
	while i < 4:
		while j < 4:
			Q = basis[i].z - basis[j].z
			Nind1 = 0
			while Nind1 < len(basis[i].alpha):
				Nind2 = 0
				while Nind2 < len(basis[j].alpha):
					p = basis[i].alpha[Nind1] + basis[j].alpha[Nind2]
					q = basis[i].alpha[Nind1] * basis[j].alpha[Nind2]/p
					P = (basis[i].alpha[Nind1]*basis[i].z + basis[j].alpha[Nind2]*basis[j].z)/p
					V_eN[i,j] += basis[i].N[Nind1]*basis[j].N[Nind2]*basis[i].c[Nind1]*basis[j].c[Nind2]*(-2.*math.pi/p)*math.exp(-1.*q*Q*Q)*(F0(p*P*P)+F0(p*(P-zH)*(P-zH)))
					Nind2 += 1
				Nind1 += 1
			j += 1
		j = 0
		i += 1
		
	return V_eN