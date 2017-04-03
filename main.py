import numpy as np
import build
import overlap
import kinetic
import eN_attract
import scf

# atomic numbers and nuclear coordinates
coords = np.matrix('0 0 0; 0 0 0.74')
V_NN = (1./coords[1,2])

### Build the basis functions
basis = build.build_basis(coords)

### overlap of the basis functions
S = overlap.overlap(basis)

### KE of an electron in each basis function
T = kinetic.kinetic(basis, S)

### electron-nucleus attraction energy
V_eN = eN_attract.attract(basis, coords)

### self-consistent field (SCF)
E, count = scf.scf_loop(T,V_eN,S,basis)
E += V_NN
print("Energy:",E,"Hartrees")
print("SCF loops:",count)