# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:07:31 2015

@author: weber
"""
import numpy as np
import matplotlib.pyplot as plt
from pressure_traces import smoothing, derivative, file_loader, filename_parse, copy, pressure_to_temperature, compress

filename = input('Filename: ')
data = file_loader(filename)
spacers, shims, Tin, pin, factor, data_date = filename_parse(filename)
pres = data[:, 1]*factor + pin*1.01325/760
time = data[:, 0]
smpr = smoothing(pres)
pc, pci = compress(smpr)
ztim = time - time[pci]
freq = np.rint(1/time[1])
dpdt = derivative(freq, pres)
smdp = smoothing(dpdt, 5)
offt = 0.001
iofig = np.argmax(smdp[(pci + offt*freq):(pci + offt*freq + 200000)])
firststage = 0.0
timeofig = time[iofig + freq*offt]*1000
timeofday = '{:02d}{:02d}'.format(data_date.hour, data_date.minute)
tempp = pres[(pci - 0.03*freq):(pci)]
temperature = pressure_to_temperature(tempp, Tin)
T_EOC = max(temperature)
copy('\t'.join(map(str, [
    timeofday, pin, Tin, pc, timeofig,
    firststage, T_EOC, spacers, shims
    ])))
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(ztim, smpr)
m = plt.get_current_fig_manager()
m.window.showMaximized()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(ztim, smpr)
ax2.plot(ztim, smdp/10000)
ax2.plot(ztim, pres)
m = plt.get_current_fig_manager()
m.window.showMaximized()
