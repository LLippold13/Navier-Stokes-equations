from parameters import *
from derivatives import *

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
###


RHS=np.ones((imax+2, jmax+2))
P*=0

while t < t_end and N<N_max:

	umax = np.amax(np.abs(U))
	vmax = np.amax(np.abs(V))

	# print('[umax, vmax] = '+str([umax,vmax]))
	if umax==0 or vmax==0:
		delt = tau * Re/(2* (1/delx/delx + 1/dely/dely))
		gamma=1
	else:
		delt = tau * min(Re/(2* (1/delx/delx + 1/dely/dely)), delx/umax, dely/vmax)
		gamma = max(umax*delt/delx,vmax*delt/dely)
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
	dUVy=duvy(U,V)
	dxxV=dxx(V)
	dyyV=dyy(V)
	dUVx=duvx(U,V)
	dqyV=dqy(V)

	# for i in np.arange(1,imax):
	# 	for j in np.arange(1,jmax):
	# 		F[i,j]=U[i,j]+delt*(1/Re(dxxU[i,j]))

	F = U + delt * (1/Re*(dxxU + dyyU) - dqxU - dUVy + GX)
	G = V + delt * (1/Re*(dxxV + dyyV) - dUVx - dqyV + GY)

	#RHS Druckgleichung


	
	for i in np.arange(1,imax+1):
		for j in np.arange(1,jmax+1):
			RHS[i,j]=1/delt*(1/delx*(F[i,j]-F[i-1,j])+1/dely*(G[i,j-G[i,j-1]]))
			if N>0:
				res[i-1,j-1] = (P[i+1,j]-2*P[i,j]+P[i-1,j])/delx/delx+(P[i,j+1]-2*P[i,j]+P[i,j-1])/dely/dely-RHS[i,j]


	#Druckberechnung
	it=0
	Pnorm = LA.norm(P,2)/(imax*jmax)
	# print('Pnorm = '+str(Pnorm))
	while it<itermax and LA.norm(res,2)/(imax*jmax) > eps*Pnorm:
		for i in np.arange(0,imax+2):
			for j in np.arange(0,jmax+2):
				if i>0 and i<imax+1 and j>0 and j<jmax+1:
					P[i,j] = (1.-omg)*P[i,j]+omg/(2*(1/delx/delx+1/dely/dely))*((P[i+1,j]+P[i-1,j])/delx/delx+(P[i,j+1]+P[i,j-1])/dely/dely-RHS[i,j])
				P[0,j]=P[1,j]
				P[i,0]=P[i,1]
				P[imax+1,j]=P[imax,j]
				P[i,jmax+1]=P[i,jmax]

				if i>0 and i<imax+1 and j>0 and j<jmax+1:
					res[i-1,j-1] = (P[i+1,j]-2*P[i,j]+P[i-1,j])/delx/delx+(P[i,j+1]-2*P[i,j]+P[i,j-1])/dely/dely-RHS[i,j]
		# print(LA.norm(res,2)/(imax*jmax) /(eps*Pnorm))
		print(it)
		it +=1

	print('it = '+str(it))
	print('[res, eps*Pnorm] = '+str([LA.norm(res,2)/(imax*jmax),(eps*Pnorm)]))

	dxP = dx(P,delx)
	dyP = dy(P,dely)
	U = F - delt * dxP
	V = G - delt * dyP

	t += delt
	N+=1
	print([t,N])
#Ende der Zeititeration
fig, ax = plt.subplots()
im = ax.imshow(P)
fig.colorbar(im)
plt.legend()
plt.show()



