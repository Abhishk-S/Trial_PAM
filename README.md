# frgd

Solver for Multiband Interacting Electron Systems via the decoupled
functional Renormalization Group. A bit unstable with only square and
hexagonal lattices implemented. Multiloop flows implemented for the vertex
but lack adequate testing. Flows at 1,2-loops have been rigorously tested
and are stable especially for moderate couplings at high temperatures.

## Description

The primary purpose of the package is the construction of the self energy and
susceptibilities of electron-election (ee) and electron-phonon (e-ph) Hamiltonians.
These quantities are crucial in constructing the ground state phase diagram of model
Hamiltonians and at the multiband level can help tackle Ab-initio Hamiltonians for
materials of interest. The decoupled fRG as implemented here is advantageous in that
it allows control over the resolution in the auxiliary frequency and momentum
channels. Such control allows for a low resolution qualitative scans of large parameter
spaces of model Hamiltonians.[^fn1]

## Installation
Older version of the package can be installed from pypi via
```console
pip install frgd
```
or the current implementation can be installed locally by running in the frgND directory.
```console
python setup.py install --user
```
## Usage

More detailed scripts can be found in the tests folder.
```python
from frgd.frg2d import frg
import numpy as np
# Initialize the parameters of the system

# Number of matsubara frequencies
nF = 20
# Temperature of the system
beTa = 32.0
# Scale at which flow is started
wMax = 100.0
# Number of sites in the system
Ns = 8
# Number of singular modes retained from SVD decomposition of bands
NB = 1
# Number of patches for expansion of the auxiliary frequency channel
NW = 4
# Number of sites for expansion of the auxiliary momentum channel
NK = 0
# Regulator for the flow: Litim, Additive, Sharp and Soft implemented
cutoffR = 'litim'
# Loop order for the fRG flow 
nL = 2
# Spatial dimension of the system
nDim = 2

# Define the current flow
currentRun = frg(nF,beTa,wMax,Ns,NW,NK,NB,cutoffR,nDim,nL)

# Momentum mode of 2D lattice
kModes = ((2*np.pi)/Ns)*np.arange(-Ns/2+1,Ns/2+1,1)
kBx = np.repeat(kModes,Ns)
kBy = np.tile(kModes,Ns)

# Number of Bands 
nB = 1

# Define single particle Hamiltonian
disMat = np.zeros((Ns*Ns,nB,nB))
disMat[:,0,0] = -2*np.cos(kBx)-2*np.cos(kBy)-4*0.15*np.cos(kBx)*np.cos(kBy)
# Define doping level (0 corresponds to half filling)
xDop = 1.0/8

# Initializes the Hamiltonian with the hopping matrix defined in disMat
currentRun.initFreeHamilt(disMat,xDop)

# Define coupling between electrons in the nB bands 
uCoup = 3
vCoup = 0.3

# Local Hubbard coupling
uIntMat = np.zeros((nB,nB))
uIntMat[0,0] = uCoup 

# Nearest neighbor density-density coupling
vIntMat = np.zeros((nB,nB))
vIntMat[0,0] = vCoup 

# Hund's coupling
jIntMat = np.zeros((nB,nB))

# Initialize electronic Interactions
currentRun.initInteractions(uIntMat,vIntMat,jIntMat)

# Optional Lines for adding a Phonon mode
# wPH = 1.0
# lEPH = 0.5
# phononDispersion = np.zeros(Ns**nDim)+wPH 
# ephCoupling = np.zeros((Ns**nDim,Ns**nDim))+ np.sqrt(wPH*lEPH/2)
# currentRun.couplePhonons(ephCoupling,phononDispersion)

# Name for output file
fileName = 'testRunOut'

# Flow with RG time l stops at lMax wCur = wMax*exp(-l)
lMax=9.0
currentRun.runFlow(fileName,lMax=lMax)

```

## Developing frgd

[^fn1]: https://arxiv.org/abs/2010.02163
