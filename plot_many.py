import matplotlib.pyplot as plt
import numpy as np
from run_many import R0_values

N_max = 7
N = len(R0_values)
N_figures = int( np.ceil( float(N) / N_max ) )

M_list = []
L_list = []
rho_list = []

i = 0
for j in range(1, N_figures + 1):
    
    plt.figure()
    plt.title("Relative convective flux for different R0 [R_sun]", fontsize=18)
    plt.grid('on')
    plt.xlabel("radius [R0]")
    plt.ylabel("ratio of F_tot")
    leg = []
    
    while i < N - (N_figures - j) * N_max:
        
        M   = []
        R   = []
        L   = []
        rho = []
        F_r = []
        F_c = []
        
        infile = open("data2/" + "r0_" + str(R0_values[i]) + ".dat", "r")
        info = infile.readline() # the first line only contains potential extra info
        for line in infile:
            line = line.split()
            M.append(float(line[1]))
            rho.append(float(line[2]))
            R.append(float(line[3]))
            L.append(float(line[5]))
            F_r.append(float(line[9]))
            F_c.append(float(line[10]))
        infile.close()
        
        M_list.append(M[-1] / M[0]) # to check how far this simulation got
        L_list.append(L[-1] / L[0])
        rho_list.append(rho[-1])
        
        R   = np.array(R) / R[0]   # [R0]
        F_r = np.array(F_r)
        F_c = np.array(F_c)
        
        F_t = F_r + F_c
        F_r /= F_t
        F_c /= F_t
        
        plt.plot(R, F_c)
        leg.append("R0=%g" % R0_values[i])
        
        i += 1

    plt.legend(leg, loc="best")

plt.show()

for i in range(N):
    print "R0 = %5.1f   M_end = %8.5f   L_end = %8.5f   rho_end = %8.1f" \
           % (R0_values[i], M_list[i], L_list[i], rho_list[i])

