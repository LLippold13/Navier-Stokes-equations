from parameters import *

import numpy as np
###

umax = np.amax(np.abs(U))
vmax = np.amax(np.abs(V))

if umax==0 and vmax==0:
	delt = tau * Re/(2* (1/delx/delx + 1/dely/dely))
else:
	delt = tau * min(Re/(2* (1/delx/delx + 1/dely/dely)), delx/umax, dely/vmax)

# Haftbedingungen
if(wW==1):
	U[0, 1:jmax+1] = 0
	V[0, 1:jmax+1] = -V[1, 1:jmax+1]
if(wO==1):
	U[imax, 1:jmax+1] = 0
	V[imax+1, 1:jmax+1] = -V[imax, 1:jmax+1]
if(wS==1):
	V[1:imax+1, 0] = 0
	U[1:imax+1, 0] = -U[1:imax+1, 1]
if(wN==1):
	V[1:imax+1, jmax] = 0
	U[1:imax+1, jmax+1] = 2*U_x -U[1:imax+1, jmax]

#Ableitungsfunktionen
def dxx(U):
	dxx=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dxx[i,j]=(U[i+1,j]-2*U[i,j]+U[i-1,j])/delx/delx
	return dxx

def dyy(U):
	dyy=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dyy[i,j]=(U[i,j+1]-2*U[i,j]+U[i,j-1])/dely/dely
	return dyy

def dqx(U):
	dqx=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dqx[i,j]=(1/delx)(np.power((U[i,j]/2+U[i+1,j]/2),2)-np.power((U[i-1,j]/2+U[i,j]/2),2))+gamma/delx*(np.abs(U[i,j]+U[i+1,j])*(U[i,j]-U[i+1,j])/4-np.abs(U[i-1,j]+U[i,j])*(U[i-1,j]-U[i,j])/4)
	return dqx

def duvy(U,V):
	duvy=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			duvy=1/dely*((V[i,j]+V[i+1,j])*(U[i,j]+U[i,j+1])/4-(V[i,j-1]+V[i+1,j-1])*(U[i,j-1]+U[i,j])/4)
			+gamma/dely*(np.abs(V[i,j]+V[i+1,j])*(U[i,j]-U[i,j+1])/4-np.abs(V[i,j-1]+V[i+1,j-1])*(U[i,j-1]-U[i,j])/4)
	return duvy

def dqy(U):
	dqy=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
		dqy[i,j]=(1/dely)(np.power((U[i,j]/2+U[i,j+1]/2),2)-np.power((U[i,j-1]/2+U[i,j]/2),2))+gamma/dely*(np.abs(U[i,j]+U[i,j+1])*(U[i,j]-U[i,j+1])/4-np.abs(U[i,j-1]+U[i,j])*(U[i,j-1]-U[i,j])/4)
	return dgy

def duvx(U,V):
	duvx=np.zeros(imax+2,jmax+2)
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			duvx=1/delx*((U[i,j]+U[i,j+1])*(V[i,j]+V[i+1,j])/4-(U[i-1,j]+U[i-1,j+1])*(V[i-1,j]+V[i,j])/4)
			+gamma/delx*(np.abs(U[i,j]+U[i,j+1])*(V[i,j]-V[i+1,j])/4-np.abs(U[i-1,j]+U[i-1,j+1])*(V[i-1,j]-V[i,j])/4)
	return duvx
 
print(delt)