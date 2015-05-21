# -*- coding: utf-8 -*-
"""
Created on Thu May 21 17:44:33 2015

@author: weber
"""
import numpy as np
import matplotlib.pyplot as plt
from pressure_traces import file_loader, filename_parse, smoothing, compress, copy
nonrfile = input('Non-reactive filename: ')
try:
    first_run
except NameError:
    first_run = True

if first_run:
    reacfile = input('Reactive filename: ')
    redata = file_loader(reacfile)
    _, _, _, rpin, rfactor, _ = filename_parse(reacfile)
    repres = redata[:, 1]*rfactor + rpin*1.01325/760
    smrepr = smoothing(repres)
    _, rpci = compress(smrepr)
    ztimere = redata[:, 0] - redata[rpci, 0]
    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    ax.plot(ztimere, smrepr)

spacers, shims, Tin, npin, nfactor, data_date = filename_parse(nonrfile)
nrdata = file_loader(nonrfile)
nrpres = nrdata[:, 1]*nfactor + npin*1.01325/760
nrsmpr = smoothing(nrpres)
maxnr = np.amax(nrsmpr)
maxnri = np.argmax(nrsmpr)
ztimenr = nrdata[:, 0] - nrdata[maxnri, 0]
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(ztimenr, nrsmpr)
timeofday = '{:02d}{:02d}'.format(data_date.hour, data_date.minute)
copy('\t'.join(map(str, [
    timeofday, npin, Tin, maxnr, 'NR', 'NR', '', spacers, shims
    ])))
