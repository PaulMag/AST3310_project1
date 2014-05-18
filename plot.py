import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 2:
    print "Usage: python plot.py name_of_datafile"
    print "Datafile is assumed to be in folder 'data2/' and have the extension '.dat'."
    sys.exit(1)

infile = open("data2/" + sys.argv[1] + ".dat", "r")

info = infile.readline() # the first line only contains potential extra info

# If you don't want to plot every step:
#resolution = 0.003
#skip = int( round( resolution / (time / n) ) )
skip = 1

dm  = []
M   = []
rho = []
R   = []
P   = []
L   = []
T   = []
eps = []
kap = []
F_r = []
F_c = []

for line in infile:
    line = line.split()

    dm.append(float(line[0]))
    M.append(float(line[1]))
    rho.append(float(line[2]))
    R.append(float(line[3]))
    P.append(float(line[4]))
    L.append(float(line[5]))
    T.append(float(line[6]))
    eps.append(float(line[7]))
    kap.append(float(line[8]))
    F_r.append(float(line[9]))
    F_c.append(float(line[10]))

infile.close()

rho0 = rho[1]
P0   = P[1]

dm  = - np.array(dm)         # [kg]
M   = np.array(M)   / M[0]   # [M0]
rho = np.array(rho) / rho[0] # [rho0]
R   = np.array(R)   / R[0]   # [R0]
P   = np.array(P)   / (1e3*P[0]) # [1e3 P0]
L   = np.array(L)   / L[0]   # [L0]
T   = np.array(T)   / 1e7    # [1ey K]
eps = np.array(eps)          # [W/kg]
kap = np.array(kap) / 1e2    # [1e3 cm^2/g]
F_r = np.array(F_r)
F_c = np.array(F_c)




# ************************************************************ #
plt.figure(); plt.hold("on")

plt.grid('on')
plt.xlabel("mass [M0]")
#plt.title("rho0 = %g kg/m^3 ,  P0 = %G Pa" % (rho0, P0), fontsize=18)

plt.plot(M[::skip], rho[::skip])
plt.plot(M[::skip], R[::skip])
plt.plot(M[::skip], P[::skip])
plt.plot(M[::skip], L[::skip])
plt.plot(M[::skip], T[::skip])
plt.plot(M[::skip], eps[::skip])
plt.plot(M[::skip], kap[::skip])
plt.legend(["rho [rho0]", \
            "R [R0]", "P [1e3 P0]", "L [L0]", "T [1e7 K]", \
            "eps [W/kg]", "kappa [1e3 cm^2/g]"], \
            loc="best")
# ************************************************************ #
plt.figure()

plt.grid('on')
plt.xlabel("mass [M0]")

plt.plot(M[::skip], dm[::skip], "b.")

plt.legend(["dm [-kg]"])
# ************************************************************ #
plt.figure()

plt.grid('on')
plt.xlabel("mass [M0]")

plt.plot(M[::skip], F_r[::skip])
plt.plot(M[::skip], F_c[::skip])

plt.legend(["F_rad", "F_con"])
# ************************************************************ #

plt.show()

