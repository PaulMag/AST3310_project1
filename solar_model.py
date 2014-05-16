# Package imports:
import os
import numpy as np
import sys
import matplotlib.pyplot as plt


# Physical constants:
c     = 2.998e8         # [m/s]
sigma = 5.67e-8         # [W/m**2/K**4]
a     = 4 * sigma / c   # [J/m**3/K**4] # constant needed to calculate T_rad
u     = 1.660538921e-27 # [kg]          # atomic mass unit
k     = 1.3806488e-23   # [J/K]
G     = 6.67384e-11     # [m**3/kg/s**2]
MeV   = 1.602176565e-13 # [J]
avogadro_inverse = 1 / 6.0221413e23 # the inverse of Avogadros number


# Ratios of different particles for star:
X     = 0.7   # hydrogen
Y_3   = 1e-10 # helium3
Y     = 0.29  # helium (total)
Z     = 0.01  # all metals
Z_7Li = 1e-13 # lithium7
Z_7Be = 1e-13 # beryllium7

mu = 1. / (2*X + 3*Y/4. + Z/2.) # average particle mass in atomic mass units


# Convection parameters:
nabla_ad = 2. / 5 # adiabatic temperature gradient for ideal gas
alpha    = 1.     # mixing length
delta    = 1.     # ideal gases
c_P      = 5/2. * k / (mu * u) # heat capacity for constant pressure


# Initial physical parameters for star:
R0   = 1.0 * 6.96e8   # [m]       # x * R_sun
L0   = 1.0 * 3.846e26 # [w]       # x * L_sun
M0   = 1.0 * 1.989e30 # [kg]      # x * M_sun
T0   = 5770           # [K]       # T_eff_sun
rho0 = 4.0e-4         # [kg/m**3]
P_rad0 = a / 3. * T0**4           # [Pa]
P_gas0 = rho0 * k * T0 / (mu * u) # [Pa] # equation of state
P0     = P_gas0 + P_rad0          # [Pa]


# Q's:
# How much energy each of the reaction give to the star when they happen.
Q_p_p     =  1.177 * MeV
Q_d_p     =  5.494 * MeV
Q_He3_He3 = 12.860 * MeV

Q_He3_He4 =  1.586 * MeV
Q_Be7_e   =  0.049 * MeV
Q_Li7_p   = 17.346 * MeV

Q_Be7_p   =  0.137 * MeV
Q_B8      =  8.367 * MeV
Q_Be8     =  2.995 * MeV

# Names which are easier to remember:
# the numbers after Q_ tells which PP chains uses this Q value
Q_123 = Q_p_p + Q_d_p # proton-deuterium reaction happens immediately when deuterium available
Q_1   = Q_He3_He3
Q_23  = Q_He3_He4
Q_2a  = Q_Be7_e
Q_2b  = Q_Li7_p
Q_3   = Q_Be7_p + Q_B8 + Q_Be8 # no reaction rates for Q_B8 and Q_Be8 available


# Read opacity file:
infile = open("opacity.txt", "r")

logR_list = infile.readline()     # list of the R's
logR_list = logR_list.split()[1:] # remove the word to the left

infile.readline() # skip empty line

logT_list = []
kappa_table = []

for line in infile: # read rest of file line by line
    line = line.split() # divide the numbers into a list
    logT_list.append(line.pop(0)) # place the first numer into the list of T's
    kappa_table.append(line) # add the rest to a nested list of kappa values
infile.close()

logT_list, logR_list, kappa_table = np.array(logT_list), np.array(logR_list), np.array(kappa_table)
logT_list   = logT_list.astype(float)
logR_list   = logR_list.astype(float)
kappa_table = kappa_table.astype(float)
kappa_x = len(logR_list)
kappa_y = len(logT_list)


# Numerical parameters:
i = 0
resolution = 5000 # how often to show progress and write to file

dm_min = - 1e15
dm_max = - 1e25
dm = dm_min

diff_min = 0.001
diff_max = 0.01


# Make placeholders for variables:
L   = L0
R   = R0
T   = T0
P   = P0
M   = M0
rho   = rho0
P_rad = P_rad0
P_gas = P_gas0

# Particle list:
class Particle:
    # class for storing simple information about each particle type

    def __init__(s, mass, ratio):
        # particle mass, particle ratio, relative particle density
        # n       = rho_rel * rho
        # rho_rel = n / rho
        if mass > 0.5:
            mass *= u # assume mass was given in [u] and convert to [kg]
    
        s.mass  = mass  # particle mass
        s.ratio = ratio # total mass ratio
        s.rho_rel = ratio / mass

H   = Particle(1.6738e-27, X)
He3 = Particle(5.0081e-27, Y_3)
He4 = Particle(6.6464e-27, Y - Y_3)
Li7 = Particle(7.01600455, Z_7Li)
Be7 = Particle(7.01692983, Z_7Be)
# Make the electron and set relative particle density to n_e = n_H + 2 * n_He:
e_  = Particle(9.10938291e-31, 0)
e_.rho_rel = H.rho_rel + 2 * He3.rho_rel + 2 * He4.rho_rel


# Functions:
def lam(i, j, T):
    # Takes two particle types and temperature as argument.
    # Checks what type of particles was given and picks the correct version of lambda.
    # Dividing my Avogadros number happens outside the function, in the loop.
    T = T / 1e9

    # PP I & II & III:
    if i == H and j == H:
        s = 4.01e-15 * T**(-2/3.) * np.exp(- 3.380 * T**(-1/3.)) \
            * ( 1 + 0.123 * T**(1/3.) + 1.09 * T**(2/3.) + 0.938 * T )
    
    # PP I:
    if i == He3 and j == He3:
        s = 6.04e10 * T**(-2/3.) * np.exp(- 12.276 * T**(-1/3.)) \
            * ( 1 + 0.034 * T**(1/3.) - 0.522 * T**(2/3.) - 0.124 * T \
               + 0.353 * T**(4/3.) + 0.213 * T**(-5/3.) )
    
    # PP II & III:
    if (i == He3 and j == He4) or (i == He4 and j == He3):
        T_star = T / (1 + 4.95e-2 * T)
        s = 5.61e6 * T_star**(5/6.) * T**(-3/2.) \
            * np.exp(- 12.826 * T_star**(-1/3.))

    # PP II:
    if (i == Be7 and j == e_) or (i == e_ and j == Be7):
        if T * 1e3 < 1: # T < 1e6
            s = 1.57e-7 / (e_.rho_rel * rho)
        else:
            s = 1.34e-10 * T**(-1/2.) * (1 - 0.537 * T**(1/3.) + 3.86 * T**(2/3.)
                + 0.0027 * T**(-1) * np.exp(2.515e-3 * T**(-1)))
    
    if (i == Li7 and j == H) or (i == H and j == Li7):
        T_star = T / (1 + 0.759 * T)
        s = 1.096e9 * T**(-2/3.) * np.exp(- 8.472 * T**(-1/3.)) \
            - 4.830e8 * T_star**(5/6.) * T**(-2/3.) \
            * np.exp(- 8.472 * T_star**(-1/3.)) \
            + 1.06e10 * T**(-3/2.) * np.exp(- 30.442 * T**(-1))

    # PP III:
    if (i == Be7 and j == H) or (i == H and j == Be7):
        s = 3.11e5 * T**(-2/3.) * np.exp(- 10.262 * T**(-1/3.)) \
            + 2.53e3 * T**(-3/2.) * np.exp(- 7.306 * T**(-1))

    #print "LAMBDA = ", s * avogadro_inverse * 1e-6 # for debugging
    return s * 1e-6 # convert to [m**3/s]


def rate(i, j, rho, T):
    # gives number of reactions per kg

    if i == j: # two particles of same kind
        delta = 1
    else: # different particles
        delta = 0

    # n_j * n_i / (rho * (1 + delta)) * lambda_ij:
    return i.rho_rel * j.rho_rel * rho / (1. + delta) * lam(i, j, T)


def kappa(T, rho):
    # Method for finding the correct kappa from the table.

    logT = np.log10(T)
    rho = rho * 0.001 # convert to [g/cm**3]
    logR = np.log10(rho / T * 1e-6)

    # Find the temperature index:
    i = 0 # this will be the vertical index in kappa_table
    for logT_l in logT_list: # loop through temperatures from small to big
        if logT_l > logT: # stop when a temp is found which is larger than T
            if i > 0:
                if logT - logT_list[i-1] < logT_l - logT: # check if the previous temp from the list was closer
                    i -= 1 # if so, take one step back
            break # no need to check the rest
        i += 1 # increase index counter and repeat (check next)
    
    # Do the same to find the R index:
    j = 0 # this will be the horizontal index in kappa_table
    for logR_l in logR_list:
        if logR_l > logR:
            if j > 0:
                if logR - logR_list[j-1] < logR_l - logR:
                    j -= 1
            break
        j += 1
    
    # If T or rho is so large that they are off the table, print a warning and
    # take one step back to avoid IndexOutOfBoundsError:
    if i >= kappa_y:
        i = kappa_y - 1
        #print "Warning: T may have blown the kappa table."
    
    if j >= kappa_x:
        j = kappa_x - 1
        #print "Warning: rho may have blown the kappa table."
    
    kappa = 10**kappa_table[i,j] # find value in table and convert back from log
    return kappa * 1000 # convert to [kg/m**3]


    def flux_tot():
        return L / (4 * np.pi * R*R)

    def flux_rad():
        return 0.75 * a * c * G * T*T*T*T * M / (kap * P * R*R)

    def flux_con():
        return flux_tot() - flux_rad()


# Data output:
if len(sys.argv) < 2: # if no output name given
    sys.argv.append("test") # default name
    print "Outfile 'data2/test.dat' was made."
if len(sys.argv) < 3: # if no info given
    sys.argv.append("") # set empty string

# Output file for results cmd-line-arg as filename:
outfile = open("data2/" + sys.argv[1] + ".dat", "w")
# If desirable, write information about current run on first line:
outfile.write(sys.argv[2] + "\n")


# Integration loop:
while True:


    # Parameters that can be found instantaneously:
    P_rad = a / 3. * T**4
    P_gas = P - P_rad
    rho   = P_gas * mu * u / (k * T)

    # The sum of the reaction rates times the energy they produce:
    eps =   rate(H  , H,   rho, T) * Q_123 \
          + rate(He3, He3, rho, T) * Q_1   \
          + rate(He3, He4, rho, T) * Q_23  \
          + rate(Be7, e_,  rho, T) * Q_2a  \
          + rate(Li7, H,   rho, T) * Q_2b  \
          + rate(Be7, H,   rho, T) * Q_3
    eps *=  avogadro_inverse # does this here to avoid doing it several times

    kap = kappa(T, rho) # find opacity

    # Sometimes print out current progress in terminal and outfile:
    if i % resolution == 0:
        print "dm  =", dm
        print "M   =", M     / M0,   "M0"
        print "rho =", rho   / rho0, "rho0"
        print "R   =", R     / R0,   "R0"
        print "P   =", P     / P0,   "P0"
        print "P_r =", P_rad / P0,   "P0"
        print "P_g =", P_gas / P0,   "P0"
        print "L   =", L     / L0,   "L0"
        print "T   =", T     / T0,   "T0"
        print "eps =", eps
        print "kap =", kap
        outfile.write("%g %f %g %g %g %g %g %g %g\n" \
                      % (dm, M, rho, R, P, L, T, eps, kap))
            # writes the current result to file for later plotting
    i += 1

    # Differential equations solved with Forward Euler:
    dR = 1. / (4. * np.pi * R*R * rho) * dm
    dP = - G * M / (4. * np.pi * R*R*R*R) * dm
    dL = eps * dm
    # dT not determined before convection check is done

    # Check for convection:
    nabla_rad = 3. * kap * L * P / \
                ( 64. * np.pi * sigma * T*T*T*T * G * M )
    #nabla_rad = (np.log(T - np.log(T)) / (np.log(P - np.log(P))
    # alternative?
    
    if nabla_rad > nabla_ad: # convective unstable => convection happens
        g = G * M / (R)     # gravity acceleration
        H_P = P / (g * rho) # pressure scale height
        # TODO: H_P different for ideal gas?
        
        U = 64 * sigma * T*T*T / (3 * kap * rho*rho * c_P) \
            * (H_P * g * delta)**0.5 # internal energy
        
        l_m = alpha * H_P # mixing length
        
        R = U / l_m
        K = 4 * R
        
        xi = np.roots([1./R, 1, K, nabla_ad-nabla_rad]) # 2 complex and 1 real
        xi = float( np.real( xi[np.imag(xi) == 0] ) ) # keeps the real root only
        
        nabla = xi*xi + K * xi + nabla_ad
        
        dT = nabla * T / P * dP

    else: # convective stable => no convection
        dT = - 3 * kap * L \
             / ( 256. * np.pi*np.pi * sigma \
             * R*R*R*R * T*T*T ) * dm

    
    # Dynamic mass step update:
    diff_largest = max( abs(dR/R), abs(dP/P), abs(dL/L), abs(dT/T) )
    
    if diff_largest > diff_max:
        if dm < dm_min: # comparison is "reverse" since dm is negative
            dm *= 0.9
    elif diff_largest < diff_min:
        if dm > dm_max:
            dm *= 1.1


    # Update parameters for next iteration:
    R += dR
    P += dP
    L += dL
    T += dT
    M += dm


    # Check if anything dropped below zero:
    # If this happens the system will stop making physical sense,
    # so the simulation should be stopped.
    # When it happens, print and save the last values of all parameters.
    if rho <= 0 or R <= 0 or P <= 0 or L <= 0 or P_rad <= 0 \
       or P_gas <= 0 or T <= 0:
        print "\nWARNING!\nSomething dropped below 0. Simulation stopped."
        print "dm  =", dm
        print "M   =", M  / M0,   "M0"
        print "rho =", rho   / rho0, "rho0"
        print "R   =", R  / R0,   "R0"
        print "P   =", P  / P0,   "P0"
        print "P_r =", P_rad / P0,   "P0"
        print "P_g =", P_gas / P0,   "P0"
        print "L   =", L  / L0,   "L0"
        print "T   =", T  / T0,   "T0"
        print "eps =", eps
        print "kap =", kap
        break

outfile.write("%g %f %g %g %g %g %g %g %g" \
              % (dm, M, rho, R, P, L, T, eps, kap)) # save last point
outfile.close()

