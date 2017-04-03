import numpy as np
import math
import eN_attract

def repulsion(b1, b2, b3, b4):
	V_ee = 0.
	PRE = 2*(math.pi**(2.5))
	
	Nind1 = 0
	Nind2 = 0
	Nind3 = 0
	Nind4 = 0
	
	while Nind1 < len(b1.alpha):	
		while Nind2 < len(b2.alpha):
			while Nind3 < len(b3.alpha):
				while Nind4 < len(b4.alpha):
					
					al = b1.alpha[Nind1]
					az = b1.z
					aN = b1.N[Nind1]
					ac = b1.c[Nind1]
					
					be = b2.alpha[Nind2]
					bz = b2.z
					bN = b2.N[Nind2]
					bc = b2.c[Nind2]
					
					ga = b3.alpha[Nind3]
					gz = b3.z
					gN = b3.N[Nind3]
					gc = b3.c[Nind3]
					
					de = b4.alpha[Nind4]
					dz = b4.z
					dN = b4.N[Nind4]
					dc = b4.c[Nind4]
					
					AB = az-bz
					GD = gz-dz
					
					p = al + be
					q = al*be/p
					Q = ((al*az) + (be*bz))/p
					normQ = aN*bN
					
					n = ga + de
					x = ga*de/n
					X = ((ga*gz) + (de*dz))/n
					normX = gN*dN
					
					QX = Q-X
					normQX = normQ*normX
					
					fint = QX*QX*p*n/(p+n)
					FTERM = eN_attract.F0(fint)
					CTERM = ac*bc*gc*dc
					expr = math.exp((q*AB*AB)+(x*GD*GD))
					POST = p*n*((p+n)**(.5))
					
					V_ee += (normQX*PRE*CTERM*FTERM)/(expr*POST)
					Nind4 += 1
				Nind4 = 0
				Nind3 += 1
			Nind3 = 0
			Nind2 += 1
		Nind2 = 0
		Nind1 += 1
	
	return V_ee