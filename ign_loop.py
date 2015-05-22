import numpy as np
import matplotlib.pyplot as plt
from pressure_traces import smoothing, derivative, file_loader, filename_parse, copy, pressure_to_temperature, compress
import os

# os.chdir('Y:\\RCM Data\\propane-dme\\phi=1.0\\75DME-25C3H8\\30-bar\\00-in-02-mm-311K')
print(os.getcwd())
pth = os.listdir('.')
result = []
for filename in pth:
    if not filename.startswith('NR_') and not filename.startswith('species') and not filename.startswith('data') and not filename.startswith('nbuoh') and os.path.isfile(filename):
        print(filename)
        data = file_loader(filename)
        spacers, shims, Tin, pin, factor, data_date = filename_parse(filename)
        pres = data[:, 1]*factor + pin*1.01325/760
        time = data[:, 0]
        smpr = smoothing(pres)
        pc, pci = compress(smpr)
        ztim = time - time[pci]
        freq = np.rint(1/time[1])
        dpdt = derivative(time, smpr)
        smdp = smoothing(dpdt, 5)
        offt = 0.001
        iofig = np.argmax(smdp[(pci + offt*freq):(pci + offt*freq + 200000)])
        firststage = 0.0
        timeofig = time[iofig + freq*offt]*1000
        timeofday = '{:02d}{:02d}'.format(data_date.hour, data_date.minute)
        tempp = pres[(pci - 0.03*freq):(pci)]
        temperature = pressure_to_temperature(tempp, Tin)
        T_EOC = max(temperature)
        result.append('\t'.join(map(str, [
            timeofday, pin, Tin, pc, timeofig,
            firststage, T_EOC, spacers, shims
            ])))
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)
        ax1.plot(ztim, smpr)

copy('\n'.join(result))
print(result)
