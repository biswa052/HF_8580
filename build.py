import math

# atomic natural orbitals (ANO) gaussian/primitive coefficients (EMSL) 
# used to construct the basis functions
alpha = [188.61445,28.276596,6.42483,1.815041,0.591063,0.212149,0.079891,0.027962]
cs1 = [0.00096385,0.00749196,0.03759541,0.14339498,0.3486363,0.43829736,0.16510661,0.02102287]
cs2 = [-0.0013119,-0.0103451,-0.0504953,-0.2073855,-0.4350885,-0.0247297,0.32252599,0.70727538]

class Orbital(object):
	
	def __init__(self,x,y,z,alpha,c):
		self.x = x
		self.y = y
		self.z = z
		self.alpha = alpha
		self.c = c
		Nind = 0
		self.N = []
		while Nind < len(self.alpha):
			(self.N).append((2*alpha[Nind]/math.pi)**(.75))
			Nind = Nind + 1
		
def build_orbital(x,y,z,alpha,c):
	orbital = Orbital(x,y,z,alpha,c)
	return orbital

def build_basis(coords):# coords.item((0,1)) = coords_1,2

	e1_1s = Orbital(coords.item((0,0)),coords.item((0,1)),coords.item((0,2)),alpha,cs1)
	e1_2s = Orbital(coords.item((0,0)),coords.item((0,1)),coords.item((0,2)),alpha,cs2)
	e2_1s = Orbital(coords.item((1,0)),coords.item((1,1)),coords.item((1,2)),alpha,cs1)
	e2_2s = Orbital(coords.item((1,0)),coords.item((1,1)),coords.item((1,2)),alpha,cs2)
	
	basis = []
	basis.append(e1_1s)
	basis.append(e1_2s)
	basis.append(e2_1s)
	basis.append(e2_2s)

	return basis