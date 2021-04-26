from matplotlib import pyplot as plt
from imitationSignal import ImSignal, Radar, Target, SClutter
from myfilter import Filter
import csv, os
import numpy as np

# Инициализация РЛС
radar1 = Radar(0.0025, 0.0003, 0.0022, 1500000, 200000000, "Rect", 4, 36, 200, 40000/(60*5))
radar2 = Radar(0.0025, 0.0003, 0.0022, 1500000, 200000000, "LFM", 4, 36, 200, 40000/(60*5))

# Инициализация целей и мешающих отражений
target1 = Target(25, 100, 102000, "gaussian", 0.1, 10, 100)
target2 = Target(30, 110, 53000, "exponential", 0.02, 20, 100)
target3 = Target(30, 89, 357800, "exponent_parabolic_1", 0.07, 30, 100)
target4 = SClutter(15, 15, 100000, 120000, "exponent_parabolic_2", 0.01, 20, 30, 10)

targetList = [target1, target2, target3, target4]

#oblako = [target4]

# Инициализация фильтра
fi1lter = Filter(0.0025, 0.1, isNormal = True)

# Инициализация объектов класса ImSignal для дальнейшего рассчёта развёрток дальности
imSignal_1 = ImSignal(radar1, targetList, fi1lter)
imSignal_2 = ImSignal(radar2, targetList, fi1lter)
#imSignal_2 = ImSignal(radar, oblako, fi1lter)

# Отчистка файла для записи развёрток дальности
with open('sweepRange.csv', 'w') as f:
	f.truncate()

# Массив, каждая строка которого, является развёрткой дальности в конкретный период повторения
allSweepRanges_1 = np.array(imSignal_1.main(500))
allSweepRanges_2 = np.array(imSignal_2.main(500))

# Построение графика азимутальной пачки в случае ЛЧМ сигнала
#-----------------------------------------------
x1 = [i for i in range(len(allSweepRanges_2))]
y1_1, y1_2, y1_3 = list(), list(), list()
allSweepRangesT = allSweepRanges_2.transpose()
for i in range(len(allSweepRanges_2)):
	y1_1.append(allSweepRangesT[1020][i])
	y1_2.append(allSweepRangesT[530][i])
	y1_3.append(allSweepRangesT[3578][i])
fig1, ax1 = plt.subplots()
ax1.set_title('Azimuth graph (LFM signal)')
ax1.step(x1, y1_1, where='post')
ax1.step(x1, y1_2, where='post')
ax1.step(x1, y1_3, where='post')
ax1.set_xlabel("Tn (repetition period number)")
ax1.set_ylabel("Signal")
ax1.legend()
ax1.grid(True)
#-----------------------------------------------
# Построение графика развёртки дальности в случае ЛЧМ сигнала
x2 = [i for i in range(len(allSweepRanges_2[0]))]
y2_1 = list(map(lambda x: x.real, allSweepRanges_2[100]))
y2_2 = list(map(lambda x: x.imag, allSweepRanges_2[100]))
y2_3 = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), allSweepRanges_2[100]))
fig2, ax2 = plt.subplots()
ax2.set_title('Range sweep graph (LFM signal)')
ax2.step(x2, y2_1, where='post')
ax2.step(x2, y2_2, where='post')
ax2.step(x2, y2_3, where='post')
ax2.set_xlabel("Range sample number")
ax2.set_ylabel("Signal")
ax2.legend()
ax2.grid(True)
#-----------------------------------------------
# Построение графика азимутальной пачки в случае простой прямоугольной модуляции сигнала
#-----------------------------------------------
x3 = [i for i in range(len(allSweepRanges_1))]
y3_1, y3_2, y3_3 = list(), list(), list()
allSweepRangesT = allSweepRanges_1.transpose()
for i in range(len(allSweepRanges_1)):
	y3_1.append(allSweepRangesT[1020][i])
	y3_2.append(allSweepRangesT[530][i])
	y3_3.append(allSweepRangesT[3578][i])
fig1, ax1 = plt.subplots()
ax3.set_title('Range sweep graph (simple rectangular modulation)')
ax3.step(x3, y3_1, where='post')
ax3.step(x3, y3_2, where='post')
ax3.step(x3, y3_3, where='post')
ax3.set_xlabel("Tn (repetition period number)")
ax3.set_ylabel("Signal")
ax3.legend()
ax3.grid(True)
#-----------------------------------------------
# Построение графика развёртки дальности в случае простой прямоугольной модуляции сигнала
#-----------------------------------------------
x4 = [i for i in range(len(allSweepRanges_1[0]))]
y4_1 = list(map(lambda x: x.real, allSweepRanges_1[100]))
y4_2 = list(map(lambda x: x.imag, allSweepRanges_1[100]))
y4_3 = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), allSweepRanges_1[100]))
fig4, ax4 = plt.subplots()
ax4.set_title('Range sweep graph (simple rectangular modulation)')
ax4.step(x4, y4_1, where='post')
ax4.step(x4, y4_2, where='post')
ax4.step(x4, y4_3, where='post')
ax4.set_xlabel("Range sample number")
ax4.set_ylabel("Signal")
ax4.legend()
ax4.grid(True)
#-----------------------------------------------
# Построение итогового графика имитатора
setrow = list()
for i in range(len(allSweepRanges_2)):
	setrow.append(list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), allSweepRanges_2[i])))
#setrow_2 = list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), imSignal_2.main(500)))
plt.matshow(setrow)
plt.colorbar()
#-----------------------------------------------
plt.show()