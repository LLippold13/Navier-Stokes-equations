import numpy as np


# geometry
xlength = 10.0
ylength = 10.0
imax = 80 	# x -> j
jmax = 80 	# y -> i
delx = xlength / (imax+1)
dely = ylength / (jmax+1)
# time iteration
t = 0.0
t_end = 1.0
delt = 1.0
tau = 0.5
N=0
N_max=3

# pressure
itermax = 100
it = 0
res=np.ones((imax, jmax))*5e8
eps = 0.01
omg = 1.7 # 0 <= omg <= 2
gamma = 1

# physics
Re = 10.0
GX = 0.0; GY = -0.0
U = np.zeros((imax+2, jmax+2)); V = np.zeros((imax+2, jmax+2)); P = np.ones((imax+2, jmax+2))
wW = 1; wO = 1; wN = 1; wS = 1
U_x = 1
F = np.zeros((imax+2, jmax+2))
G = np.zeros((imax+2, jmax+2))
