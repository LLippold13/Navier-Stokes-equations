from parameters import *

import numpy as np
###

#Ableitungsfunktionen
def dxx(U):
	dxx=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dxx[i,j]=(U[i+1,j]-2*U[i,j]+U[i-1,j])/delx/delx
	return dxx

def dyy(U):
	dyy=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dyy[i,j]=(U[i,j+1]-2*U[i,j]+U[i,j-1])/dely/dely
	return dyy

def dqx(U):
	dqx=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dqx[i,j]=(1/delx)*(np.power((U[i,j]/2+U[i+1,j]/2),2)-np.power((U[i-1,j]/2+U[i,j]/2),2))+gamma/delx*(np.abs(U[i,j]+U[i+1,j])*(U[i,j]-U[i+1,j])/4-np.abs(U[i-1,j]+U[i,j])*(U[i-1,j]-U[i,j])/4)
	return dqx

def duvy(U,V):
	duvy=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			duvy=1/dely*((V[i,j]+V[i+1,j])*(U[i,j]+U[i,j+1])/4-(V[i,j-1]+V[i+1,j-1])*(U[i,j-1]+U[i,j])/4)
			+gamma/dely*(np.abs(V[i,j]+V[i+1,j])*(U[i,j]-U[i,j+1])/4-np.abs(V[i,j-1]+V[i+1,j-1])*(U[i,j-1]-U[i,j])/4)
	return duvy

def dqy(U):
	dqy=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dqy[i,j]=(1/dely)*(np.power((U[i,j]/2+U[i,j+1]/2),2)-np.power((U[i,j-1]/2+U[i,j]/2),2))+gamma/dely*(np.abs(U[i,j]+U[i,j+1])*(U[i,j]-U[i,j+1])/4-np.abs(U[i,j-1]+U[i,j])*(U[i,j-1]-U[i,j])/4)
	return dqy

def duvx(U,V):
	duvx=np.zeros((imax+2,jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			duvx=1/delx*((U[i,j]+U[i,j+1])*(V[i,j]+V[i+1,j])/4-(U[i-1,j]+U[i-1,j+1])*(V[i-1,j]+V[i,j])/4)
			+gamma/delx*(np.abs(U[i,j]+U[i,j+1])*(V[i,j]-V[i+1,j])/4-np.abs(U[i-1,j]+U[i-1,j+1])*(V[i-1,j]-V[i,j])/4)
	return duvx

def dx(P, delx):
	dx = np.zeros((imax+2, jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dx = (P[i+1,j] - P[i,j])/delx
	return dx

def dy(P, dely):
	dy = np.zeros((imax+2, jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			dy = (P[i,j+1] - P[i,j])/dely
	return dy



while t < t_end:

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



	 # F und G berechnen
	dxxU=dxx(U)
	dyyU=dyy(U)
	dqxU=dqx(U)
	duvy=duvy(U,V)
	dxxV=dxx(V)
	dyyV=dyy(V)
	duvx=duvx(U,V)
	dqyV=dqy(V)

	# for i in np.arange(1,imax):
	# 	for j in np.arange(1,jmax):
	# 		F[i,j]=U[i,j]+delt*(1/Re(dxxU[i,j]))

	F = U + delt * (1/Re*(dxxU + dyyU) - dqxU - duvy + GX)
	G = V + delt * (1/Re*(dxxV + dyyV) - duvx - dqyV + GY)

	#RHS Druckgleichung


	RHS=np.zeros((imax+2, jmax+2))
	for i in np.arange(1,imax):
		for j in np.arange(1,jmax):
			RHS[i,j]=1/delt*(1/delx*(F[i,j]-F[i-1,j])+1/dely*(G[i,j-G[i,j-1]]))

	#Druckberechnung
	it=0
	while it<=itermax and r < eps*np.amax(res):
		for i in np.arange(1,imax):
			for j in np.arange(1,jmax):
				P[i,j] = (1-omg)*P[i,j]+omg/(2*(1/delx/delx+1/dely/dely))*((P[i+1,j]+P[i-1,j])/delx/delx+(P[i,j+1]+P[i-1,j-1])/dely/dely-RHS[i,j])
				res[i,j] = (P[i+1,j]-2*P[i,j]+P[i-1,j])/delx/delx+(P[i,j+1]-2*P[i,j]+P[i,j-1])/dely/dely-RHS[i,j]
				P[0,j]=P[1,j]
				P[i,0]=P[i,1]
				P[imax+1,j]=P[imax,j]
				P[i,jmax+1]=P[i,jmax]

	dxP = dx(P,delx)
	dyP = dy(P,dely)
	U = F - delt * dxP
	V = G - delt * dyP

	t += delt
#Ende der Zeititeration


