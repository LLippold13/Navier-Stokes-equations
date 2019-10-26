import numpy as np


# geometry
xlength = 1.0
ylength = 1.0
imax = 10 	# x -> j
jmax = 10 	# y -> i
delx = xlength / imax
dely = ylength / jmax

# time iteration
t = 0.0
t_end = 1.0
delt = 0.1
tau = 0.5

# pressure
itermax = 10
it = 0
res = 1e10+0.0
eps = 1.0
omg = 1.7
gamma = 1.0

# physics
Re = 1.0
GX = 0.0; GY = -10.0
U = np.zeros((imax+2, jmax+2)); V = np.zeros((imax+2, jmax+2)); P = np.zeros((imax+2, jmax+2))
wW = 1; wO = 1; wN = 1; wS = 1
U_x = 1
F = np.zeros((imax+2, jmax+2))
G = np.zeros((imax+2, jmax+2))
