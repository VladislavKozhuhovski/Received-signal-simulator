from matplotlib import pyplot as plt
from imitationSignal import ImSignal, Radar, Target, SClutter
from myfilter import Filter
import csv, os
import numpy as np

radar = Radar(0.0025, 0.0003, 0.0022, 1500000, 200000000, "Rect", 4, 36)

target1 = Target(25, 100, 102000, "gaussian", 0.1, 10)
target2 = Target(30, 110, 53000, "exponential", 0.02, 20)
target3 = Target(22, 89, 357800, "exponent_parabolic_1", 0.07, 30)
target4 = SClutter(40, 15, 100000, 120000, "exponent_parabolic_2", 0.01, 20, 30)

targetList = [target1, target2, target3, target4]

oblako = [target4]

fi1lter = Filter(0.0025, 0.1, isNormal = True)

#imSignal = ImSignal(radar, targetList, fi1lter)

imSignal_2 = ImSignal(radar, oblako, fi1lter)

with open('sweepRange.csv', 'w') as f:
	f.truncate()

#setrow = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), imSignal.main(500)))

setrow_2 = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), imSignal_2.main(500)))

plt.plot(setrow_2[0])
#plt.matshow(setrow)
#plt.colorbar()
plt.show()