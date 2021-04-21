from matplotlib import pyplot as plt
from imitationSignal import ImSignal, Radar, Target, SClutter
from myfilter import Filter
import csv, os
import numpy as np

# Инициализация РЛС
radar = Radar(0.0025, 0.0003, 0.0022, 1500000, 200000000, "Rect", 4, 36)

# Инициализация целей и мешающих отражений
target1 = Target(25, 100, 102000, "gaussian", 0.1, 10)
target2 = Target(30, 110, 53000, "exponential", 0.02, 20)
target3 = Target(40, 89, 357800, "exponent_parabolic_1", 0.07, 30)
target4 = SClutter(15, 15, 100000, 120000, "exponent_parabolic_2", 0.01, 20, 30)

targetList = [target1, target2, target3, target4]

#oblako = [target4]

# Инициализация фильтра
fi1lter = Filter(0.0025, 0.1, isNormal = True)

# Инициализация объектра класса ImSignal для дальнейшего рассчёта развёрток дальности
imSignal = ImSignal(radar, targetList, fi1lter)

#imSignal_2 = ImSignal(radar, oblako, fi1lter)

# Отчистка файла для записи развёрток дальности
with open('sweepRange.csv', 'w') as f:
	f.truncate()

# Массив, каждая строка которого, является развёрткой дальности в конкретный период повторения
setrow = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), imSignal.main(500)))

#setrow_2 = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), imSignal_2.main(500)))

#plt.plot(setrow_2[0])

# Построение графика
plt.matshow(setrow)
plt.colorbar()
plt.show()