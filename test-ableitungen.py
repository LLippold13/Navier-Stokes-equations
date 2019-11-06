from parameters import *

import matplotlib.pyplot as plt


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
			duvy[i,j]=1/dely*((V[i,j]+V[i+1,j])*(U[i,j]+U[i,j+1])/4-(V[i,j-1]+V[i+1,j-1])*(U[i,j-1]+U[i,j])/4)
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
			duvx[i,j]=1/delx*((U[i,j]+U[i,j+1])*(V[i,j]+V[i+1,j])/4-(U[i-1,j]+U[i-1,j+1])*(V[i-1,j]+V[i,j])/4)
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



T=np.zeros((imax+2,jmax+2))
for i in range(imax+2):
	for j in range(jmax+2):
		U[i,j] = i*delx
		V[i,j] = i*delx
		T[i,j] = 2*i*delx
		# U[i,j] = delx*i*delx*i + dely*j*dely*j
		# T[i,j] = 2
		# U[i,j]=np.sin(i*delx)*np.cos(j*dely)
		# T[i,j]=2*np.cos(i*delx)*np.sin(i*delx)*np.cos(j*dely)*np.cos(j*dely)

dUV = duvx(U,V)
# print((dU[50,50]-T[50,50])/dU[50,50])
fig, ax = plt.subplots()
im = ax.imshow((dUV-T)/dUV)
fig.colorbar(im)
plt.legend()
plt.show()