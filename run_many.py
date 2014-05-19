import sys, os

R0_values = [0.4, 0.7, 1, 1.5, 2.5, 4, 7, \
             10, 15, 20, 25, 30, 35, 40, \
             45, 50, 55, 60, 65, 70, 75]

if __name__ == "__main__":

    i = 0

    for R0 in R0_values:
        print "%d / %d: Working on %g" % (i, len(R0_values), R0)
        os.system("python solar_model.py %g r0_%g" % (R0, R0))
        i += 1

    print "%d / %d: FINISHED!" % (i, len(R0_values))

