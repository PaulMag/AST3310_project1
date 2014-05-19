import matplotlib.pyplot as plt
import numpy as np
from run_many import R0_values

N_max = 7
N = len(R0_values)
N_figures = int( np.ceil( float(N) / N_max ) )

M_list = []

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
        F_r = []
        F_c = []
        
        infile = open("data2/" + "r0_" + str(R0_values[i]) + ".dat", "r")
        info = infile.readline() # the first line only contains potential extra info
        for line in infile:
            line = line.split()
            M.append(float(line[1]))
            R.append(float(line[3]))
            F_r.append(float(line[9]))
            F_c.append(float(line[10]))
        infile.close()
        
        M_list.append(M[-1] / M[0]) # to check how far this simulation got
        
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

M_best = min(M_list)
for i in range(N):
    print "R0 = %5.1f   M_last = %7.5f" % (R0_values[i], M_list[i])
    
    if M_list[i] == M_best:
        i_best = i

print "R0 = %g gave the best result." % R0_values[i_best]

