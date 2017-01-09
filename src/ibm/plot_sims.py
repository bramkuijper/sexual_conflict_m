#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *

# set some basic parameters
params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# open the file and get parameters
f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^type.*",line) != None:
        parline = idx - 1;
        break;

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    histdat = pd.read_csv(filename, sep=";")

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 4

if "meanoff_p" in histdat.columns.values:
    num_rows = 5

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(
        histdat["generation"],
        histdat["meanoff"],
        'darkgreen',
        linewidth=1
        )

plt.plot(
        histdat["generation"],
        histdat["meansen"],
        'red',
        linewidth=1
        )

plt.plot(
        histdat["generation"],
        histdat["meanthr"],
        'blue',
        linewidth=1
        )
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Mean phenotypes')

plt.legend((r"$\bar{p}$",r"$\bar{s}$",r"$\bar{t}$"))

# add second subplot
plt.subplot(num_rows,1,2)
plt.plot(
        histdat["generation"],
        histdat["varoff"],
        'darkgreen',
        linewidth=1
        )

plt.plot(
        histdat["generation"],
        histdat["varsen"],
        'red',
        linewidth=1
        )

plt.plot(
        histdat["generation"],
        histdat["varthr"],
        'blue',
        linewidth=1
        )

plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Variances')

plt.legend((r"$\sigma_{p}^{2}$",r"$\sigma_{s}^{2}$",r"$\sigma_{t}^{2}$"))

# add third subplot
plt.subplot(num_rows,1,3)
plt.plot(
        histdat["generation"],
        histdat["Nm"],
        'blue',
        linewidth=1)

plt.plot(
        histdat["generation"],
        histdat["Nf"],
        'red',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Population size')
plt.legend((r"$N_{\mathrm{m}}$",r"$N_{\mathrm{f}}$"))

# add fourth subplot
plt.subplot(num_rows,1,4)
plt.plot(
        histdat["generation"],
        histdat["Nmsurv"],
        'blue',
        linewidth=1)

plt.plot(
        histdat["generation"],
        histdat["Nfsurv"],
        '#ff004d',
        linewidth=1)
plt.plot(
        histdat["generation"],
        histdat["Nfunmated"],
        'darkgreen',
        linewidth=1)
plt.ylabel(r'Survivors')
plt.legend((
                r'$N_{\mathrm{m}}^{*}$',
                r'$N_{\mathrm{f}}^{*}$',
                r'$N_{\mathrm{f,unmated}}$'
                ))

# add subplot maternal and paternal effects
if "meanoff_p" in histdat.columns.values:
    plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')

    plt.subplot(num_rows,1,5)

    plt.plot(
            histdat["generation"],
            histdat["meanoff_p"],
            'darkgreen',
            linewidth=1)

    plt.plot(
            histdat["generation"],
            histdat["meanthr_m"],
            'blue',
            linewidth=1)

    plt.plot(
            histdat["generation"],
            histdat["meansen_m"],
            'red',
            linewidth=1)
    plt.ylabel(r'Transgen')

    plt.legend((
                    r'$\bar{\mathcal{P}}_{p}$',
                    r'$\bar{m}_{t}$',
                    r'$\bar{m}_{s}$',
                    ))

plt.xlabel(r'Generations')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
