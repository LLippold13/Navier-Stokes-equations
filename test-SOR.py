from parameters import *
from derivatives import *

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

###



for i in range(imax+2):
	for j in range(jmax+2):
		# P[i,j] = 2 - (delx*i + dely*j)
		# P[i,j]=np.sin(i*delx)*np.cos(j*dely)
		# P[i,j] = np.sin((delx*i + dely*j))
		P[i,j] = delx*i + dely*j

while t < t_end and N<N_max:

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
	dUVy=duvy(U,V)
	dxxV=dxx(V)
	dyyV=dyy(V)
	dUVx=duvx(U,V)
	dqyV=dqy(V)



	F = U + delt * (1/Re*(dxxU + dyyU) - dqxU - dUVy + GX)
	G = V + delt * (1/Re*(dxxV + dyyV) - dUVx - dqyV + GY)

	#RHS Druckgleichung


	RHS=np.ones((imax+2, jmax+2))
	RHS *= 4

	for i in np.arange(0,imax+2):
		for j in np.arange(0,imax+2):
			# P[i,j] = np.sin(i*delx) + np.cos(j*dely) +i*delx*i*delx
			# P[i,j] = i*delx*i*delx + j*dely*j*dely*2
			P[i,j] = 1

	#Druckberechnung
	it=0
	Pnorm = LA.norm(P,2)/(imax*jmax)
	print('Pnorm = '+str(Pnorm))
	fig, ax = plt.subplots()

	for i in np.arange(1,imax+1):
		for j in np.arange(1,jmax+1):
			res[i-1,j-1] = (P[i+1,j]-2*P[i,j]+P[i-1,j])/delx/delx+(P[i,j+1]-2*P[i,j]+P[i,j-1])/dely/dely-RHS[i,j]

	print('res = '+str(LA.norm(res,2)/(imax*jmax)))


	while it<=itermax and LA.norm(res,2)/(imax*jmax) >= eps*Pnorm:
		for i in np.arange(1,imax+1):
			P[i,jmax+1]=i*delx*i*delx + ylength*ylength
			P[i,0]=i*delx*i*delx
			for j in np.arange(1,jmax+1):
				P[0,j]=j*dely*j*dely
				P[imax+1,j]= (imax+1)*delx*(imax+1)*delx + j*dely*j*dely

				# P[0,j]=P[1,j]
				# P[i,0]=P[i,1]
				# P[imax+1,j]=P[imax,j]
				# P[i,jmax+1]=P[i,jmax]
		for i in np.arange(1,imax+1):
			for j in np.arange(1,jmax+1):
				P[i,j] = (1.-omg)*P[i,j]+omg/(2.*(1./(delx*delx)+1./(dely*dely))) * ((P[i+1,j]+P[i-1,j])/(delx*delx)+(P[i,j+1]+P[i,j-1])/(dely*dely)-RHS[i,j])
		for i in np.arange(1,imax+1):
			for j in np.arange(1,jmax+1):
				res[i-1,j-1] = (P[i+1,j]-2*P[i,j]+P[i-1,j])/delx/delx+(P[i,j+1]-2*P[i,j]+P[i,j-1])/dely/dely-RHS[i,j]
		# plt.pause(0.001)
		print('res/eps = '+str(LA.norm(res,2)/(imax*jmax) /(Pnorm*eps)))
		print(it)
		it +=1
		


	dxP = dx(P,delx)
	dyP = dy(P,dely)
	U = F - delt * dxP
	V = G - delt * dyP

	t += delt
	N+=1
	print([t,N])
#Ende der Zeititeration



im = ax.imshow(P)
plt.figure()
plt.plot(np.arange(0,imax+2), np.transpose(P[:,20]))
fig.colorbar(im)
plt.legend()



plt.show()




