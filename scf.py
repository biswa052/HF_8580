import numpy as np
import ee_repulsion4
import copy

def GfromP(P,basis):
	G = np.zeros((4,4))
	x = 0
	y = 0
	i = 0
	j = 0
	while x<4:
		while y<4:
			while i<4:
				while j<4:
					gxyij = ee_repulsion4.repulsion(basis[x],basis[y],basis[i],basis[j])
					G[x,y] += P[i,j]*gxyij
					j+=1
				j=0
				i+=1
			y+=1
			i=0
		y=0
		x+=1
	return G

def XCfromP(P,basis):
	XC = np.zeros((4,4))
	while x<4:
		while y<4:
			XC[x,y] += -1.*1.25*(10.**-4.)*vxxy(P,basis,x,y) #
			#XC[x,y] += vcxy(P,basis,x,y) contribution from ec
			#XC[x,y] += vcpxy(P,basis,x,y) contribiutin from ec'
			y+=1
		y=0
		x+=1
	return XC

def PfromC(C,eps):
	P = np.zeros((4,4))
	i = 0
	j = 0
	a = min(range(len(eps)), key=eps.__getitem__)
	while i < 4:
		while j < 4:
			P[i,j] = 2.*(C[i,a]*np.conjugate(C[j,a]))
			j+=1
		j=0
		i+=1
	return P
	
def calc_E(P,H_core,F):
	i = 0
	j = 0
	E = 0.
	while i < 4:
		while j < 4:
			E += P[j,i]*(H_core[i,j] + F[i,j])
			j += 1
		j = 0
		i += 1
	E *= .5
	return E

def scf_loop(T,V_eN,S,basis):
	H_core = T + V_eN
	
	s, U = np.linalg.eigh(S)
	Ut = np.transpose(U)
	shalf = np.zeros((4,4))
	sint = 0
	while sint < 4:
		shalf[sint,sint] = s[sint]**(-.5)
		sint += 1
	X = np.dot(U,shalf)
	X = np.dot(X,Ut)
	Xt = np.transpose(X)
	
	P = np.zeros((4,4)) # guess for density matrix 
	F = H_core
	Eold = calc_E(P,H_core,F)
	
	Fp = np.dot(Xt,F)
	Fp = np.dot(Fp,X)
	eps, Cp = np.linalg.eigh(Fp)
	C = np.dot(X,Cp)
	P = PfromC(C,eps)
			
	check = 1
	count = 0
	while(check==1 or (abs(Enew-Eold) > (10.**(-8)))):
		if(check==0):
			Eold = Enew
			
		G = GfromP(P,basis)
		F = H_core + G
		Enew = calc_E(P,H_core,F)
		
		Fp = np.dot(Xt,F)
		Fp = np.dot(Fp,X)
		
		eps, Cp = np.linalg.eig(Fp)
		C = np.dot(X,Cp)
		P = PfromC(C,eps)
				
		check = 0
		count += 1
	print(P)
	return Enew, count