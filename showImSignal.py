from matplotlib import pyplot as plt
from imitationSignal import ImSignal, Radar, Target, SClutter
from myfilter import Filter
import csv, os

radar = Radar(0.0025, 0.0003, 0.0022, 1500000, 200000000, "Rect", 4, 36)

target1 = Target(25, 100, 102000, "gaussian", 0.1, 40)
target2 = Target(30, 110, 53000, "exponential", 0.02, 55)
target3 = Target(22, 89, 357400, "exponent_parabolic_1", 0.07, 222)
target4 = SClutter(15, 15, 240900, 357400, "exponent_parabolic_2", 0.01, 202, 222)

targetList = [target1, target2, target3, target4]

fi1lter = Filter(0.0025, 0.1, isNormal = True)

imSignal = ImSignal(radar, targetList, fi1lter)

with open('sweepRange.csv', 'w') as f:
	f.truncate()

with open('sweepRange.csv', 'r') as f:
	imSignal.main(10)
	sweep_reader = csv.reader(f)
	for row in sweep_reader:
		if row:
			x = [i for i in range(len(row))]
			plt.scatter(x, row)
			plt.show()
			break