
import numpy as np 
from frgd.frg2d import frg 

nF=24
wMax=100.0
NB=1
nDim=3
cutoffR='litim'

nL=2
beTa=40
uU=6
V0=0.66
V1=0.84
V2=0
xDopC=0
xDopF=0.1
NW=4
NK=0
Ns=8
outputName='frg2DPAT0N8U3v00v10'


currentRun=frg(nF,beTa,wMax,Ns,NW,NK,NB,cutoffR,nDim,nL)
nB=1
kX=((2*np.pi)/float(Ns))*np.arange(-Ns/2+1,Ns/2+1,1)
kBx=np.tile(kX[:,None,None],(1,Ns,Ns))
kBy=np.tile(kX[None,:,None],(Ns,1,Ns))
kBz=np.tile(kX[None,None,:],(Ns,Ns,1))

disMat=np.zeros((Ns,Ns,Ns,nB,nB))
disMat[:,:,:,0,0]=-2*np.cos(kBx)-2*np.cos(kBy)-2*np.cos(kBz)
disMat=np.reshape(disMat,(Ns*Ns*Ns,nB,nB))

vCF=np.zeros((Ns,Ns,Ns,nB,nB))
vCF[:,:,:,0,0]=V0+2*V1*(np.cos(kBx)+np.cos(kBy)+np.cos(kBz))+4*V2*np.cos(kBx)*np.cos(kBy)*np.cos(kBz)
vCF=np.reshape(vCF,(Ns*Ns*Ns,nB,nB))

uIntMat=np.zeros((nB,nB))
uIntMat[0,0]=uU

currentRun.initFreePeriodicHamilt(disMat,xDopC,vCF,xDopF)
currentRun.initInteractions(uIntMat)
fileName=outputName
currentRun.runFlow(fileName)
