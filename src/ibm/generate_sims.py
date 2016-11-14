#!/usr/bin/env python

import numpy as np

exename = "./xsexconfl"


a = [0.01, 0.1, 0.5, 2, 5 ]
co = [ 0.001, 0.01, 0.1, 0.2 ]
ct = [ 0.001, 0.01, 0.1, 0.2 ]
cs = [ 0.001, 0.01, 0.1, 0.2 ]
theta_psi = [ 0.5, 1.0, 2.0, 0.1 ]

ctr = 0

for a_i in a:
    for co_i in co:
        for ct_i in ct:
            for cs_i in cs:
                for theta_psi_i in theta_psi:

                    ctr += 1
                    print("echo " + str(ctr))
                    print(exename + " " + str(a_i)
                            + " " + str(co_i) + " " + str(ct_i)
                            + " " + str(cs_i)
                            + " 0.01 0.01 0.01 0.02 0.02 0.02"
                            + " " + str(theta_psi_i))
#	a = atof(argv[1]);
#	co = atof(argv[2]);
#	ct = atof(argv[3]);
#	cs = atof(argv[4]);
#	mu_off = atof(argv[5]);
#	mu_thr = atof(argv[6]);
#	mu_sen = atof(argv[7]);
#	sdmu_off = atof(argv[8]);
#	sdmu_thr = atof(argv[9]);
#	sdmu_sen = atof(argv[10]);
#    theta_psi = atof(argv[11]);
