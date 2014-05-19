import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 2:
    print "Usage: python plot.py name_of_datafile"
    print "Datafile is assumed to be in folder 'data2/' and have the extension '.dat'."
    sys.exit(1)

infile = open("data2/" + sys.argv[1] + ".dat", "r")

info = infile.readline() # the first line only contains potential extra info

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
eprod_123 = []
eprod_1   = []
eprod_2   = []
eprod_3   = []

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
    eprod_123.append(float(line[11]))
    eprod_1.append(float(line[12]))
    eprod_2.append(float(line[13]))
    eprod_3.append(float(line[14]))

infile.close()


dm  = - np.array(dm) # SI-units
M   = np.array(M)
rho = np.array(rho)
R   = np.array(R)
P   = np.array(P)
L   = np.array(L)
T   = np.array(T)
eps = np.array(eps)
kap = np.array(kap)

F_r = np.array(F_r)
F_c = np.array(F_c)
F_t = F_r + F_c

eprod_123 = np.array(eprod_123)
eprod_1   = np.array(eprod_1)
eprod_2   = np.array(eprod_2)
eprod_3   = np.array(eprod_3)

R_sun = 6.96e8
R_info = "R0 = %4.1f R_sun" % (R[0] / R_sun)
R /= R[0] # [R0]

# ************************************************************ #
plt.figure(); plt.hold("on")
plt.title("Parameters,  " + R_info, fontsize=18)

plt.grid('on')
plt.xlabel("radius [R0]")

plt.plot(R, M / M[0])
plt.plot(R, np.log10(P / (1e3 * P[0])))
plt.plot(R, L / L[0])
plt.plot(R, np.log10(T / 1e3))

plt.legend(["M [M0]", "P [log10(1e3 P0)]", "L [L0]", "T [log10(1e3 K])"], \
            loc="best")
# ************************************************************ #
plt.figure()
plt.title("Mass step distribution,  " + R_info, fontsize=18)

plt.grid('on')
plt.xlabel("radius [R0]")

plt.plot(R, np.log10(dm), "b.")

plt.legend(["dm [log10(-kg)]"], loc="best")
# ************************************************************ #
plt.figure()
plt.title("Flux relations,  " + R_info, fontsize=18)

plt.grid('on')
plt.xlabel("radius [R0]")
plt.ylabel("ratio of F_tot")

plt.plot(R, F_r / F_t)
plt.plot(R, F_c / F_t)

plt.legend(["F_rad", "F_con"], loc="best")
# ************************************************************ #
plt.figure()
plt.title("Energy production of PP-chains,  " + R_info, fontsize=18)

plt.grid('on')
plt.xlabel("radius [R0]")
plt.ylabel("log10(N_a * W / kg)")
plt.axis([0, 1, -50, 25]) # this may be adjusted for each situation

plt.plot(R, np.log10(eprod_123))
plt.plot(R, np.log10(eprod_1))
plt.plot(R, np.log10(eprod_2))
plt.plot(R, np.log10(eprod_3))

plt.legend(["pp + pd", "PP1", "PP2", "PP3"], loc="best")
# ************************************************************ #

plt.show()

