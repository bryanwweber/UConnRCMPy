# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:13:12 2015

@author: weber
"""
import numpy as np
import matplotlib.pyplot as plt
from pressure_traces import smoothing, compress, pressure_to_volume, pressure_to_temperature, volume_to_pressure, filename_parse

nonrfile = 'NR_00_in_00_mm_311K-1115t-40x-15-May-15-1425.txt'
reacfile = '00_in_00_mm_311K-1118t-40x-15-May-15-1301.txt'
nonrdata = np.fromfile(nonrfile, sep=' ')
nonrdata = np.reshape(nonrdata, (len(nonrdata)/2, 2))
_, _, _, npin, nfactor, _ = filename_parse(nonrfile)
reacdata = np.fromfile(reacfile, sep=' ')
reacdata = np.reshape(reacdata, (len(reacdata)/2, 2))
_, _, reacTin, rpin, rfactor, _ = filename_parse(reacfile)
reacoffs = 0
nonroffs = 0
comptime = 32
nonrend = 400
reacend = 80
reactime = reacdata[:, 0]
reacpres = reacdata[:, 1]*rfactor + rpin*1.01325/760
nonrtime = nonrdata[:, 0]
nonrpres = nonrdata[:, 1]*nfactor + npin*1.01325/760

reacsmpr = smoothing(reacpres)
nonrsmpr = smoothing(nonrpres)
reacpc, reacpci = compress(reacsmpr)
nonrpci = np.argmax(nonrsmpr)
reacztim = reactime - reactime[reacpci]
reacfreq = np.rint(1/reactime[1])
nonrfreq = np.rint(1/nonrtime[1])
nonrendidx = nonrend/1000*nonrfreq
reacendidx = reacend/1000*reacfreq
startpoint = comptime/1000*reacfreq
stroke_pressure = reacsmpr[(reacpci - startpoint):(reacpci + reacoffs)]
post_pressure = nonrsmpr[(nonrpci + nonroffs):(nonrpci +
                         nonrendidx + nonroffs - reacoffs)]
print_pressure = reacsmpr[(reacpci - startpoint):(reacpci + reacendidx)]
time = np.arange(-comptime/1000, nonrend/1000, 1/reacfreq)
stroke_volume = pressure_to_volume(stroke_pressure, reacTin)
stroke_temperature = pressure_to_temperature(stroke_pressure, reacTin)
compressed_volume = stroke_volume[-1]
compressed_temperature = stroke_temperature[-1]
post_volume = pressure_to_volume(post_pressure,
                                 compressed_temperature,
                                 compressed_volume)
volume = np.concatenate((stroke_volume, post_volume))
pressure = volume_to_pressure(volume, stroke_pressure[0]*1E5, reacTin)

fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
ax.plot(reacztim, reacsmpr)
ax.plot(time[0:len(print_pressure)], print_pressure)
ax.plot(time, pressure)
m = plt.get_current_fig_manager()
m.window.showMaximized()

volout = np.vstack((time[::5] + comptime/1000, volume[::5])).transpose()
presout = np.vstack((time[:len(print_pressure):5] + comptime/1000, print_pressure[::5])).transpose()
np.savetxt('volume.csv', volout, delimiter=',')
np.savetxt('Tc__P0__T0_{}K_pressure.txt'.format(reacTin), presout, delimiter='\t')
