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
target2 = Target(30, 110, 53000, "exponential", 0.02, 10, 100)
target3 = Target(30, 89, 357800, "exponent_parabolic_1", 0.07, 10, 100)
target4 = SClutter(15, 15, 100000, 120000, "exponent_parabolic_2", 0.01, 20, 30, 10)

targetList = [target1, target2, target3]

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

# Построение графиков азимутальных пачек в случае ЛЧМ сигнала
#-----------------------------------------------
x1 = [i for i in range(len(allSweepRanges_2))]
y1_1, y1_2, y1_3 = list(), list(), list()
allSweepRangesT = allSweepRanges_2.transpose()
for i in range(len(allSweepRanges_2)):
	y1_1.append(allSweepRangesT[1020][i])
	y1_2.append(allSweepRangesT[530][i])
	y1_3.append(allSweepRangesT[3578][i])
# Первая цель
fig1_1, ax1_1 = plt.subplots()
ax1_1.set_title('Azimuth graph (LFM signal), target 1')
ax1_1.step(x1, list(map(lambda x: x.real, y1_1)), where='post', label = 'Re')
ax1_1.step(x1, list(map(lambda x: x.imag, y1_1)), where='post', label = 'Im')
ax1_1.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y1_1)), where='post', label = 'Mod')
ax1_1.set_xlabel("Tn (repetition period number)")
ax1_1.set_ylabel("Signal")
ax1_1.legend()
ax1_1.grid(True)
# Спектр азимутальной пачки
fig1_1_2, ax1_1_2 = plt.subplots()
ax1_1_2.set_title('Spectrum of Azimuth graph (LFM signal), target 1')
ax1_1_2.step(x1, abs(np.fft.fft(y1_1)))
ax1_1_2.legend()
ax1_1_2.grid(True)
# Вторая цель
fig1_2, ax1_2 = plt.subplots()
ax1_2.set_title('Azimuth graph (LFM signal), target 2')
ax1_2.step(x1, list(map(lambda x: x.real, y1_2)), where='post', label = 'Re')
ax1_2.step(x1, list(map(lambda x: x.imag, y1_2)), where='post', label = 'Im')
ax1_2.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y1_2)), where='post', label = 'Mod')
ax1_2.set_xlabel("Tn (repetition period number)")
ax1_2.set_ylabel("Signal")
ax1_2.legend()
ax1_2.grid(True)
# Спектр азимутальной пачки
fig1_2_2, ax1_2_2 = plt.subplots()
ax1_2_2.set_title('Spectrum of Azimuth graph (LFM signal), target 2')
ax1_2_2.step(x1, abs(np.fft.fft(y1_2)))
ax1_2_2.legend()
ax1_2_2.grid(True)
# Третья цель
fig1_3, ax1_3 = plt.subplots()
ax1_3.set_title('Azimuth graph (LFM signal), target 3')
ax1_3.step(x1, list(map(lambda x: x.real, y1_3)), where='post', label = 'Re')
ax1_3.step(x1, list(map(lambda x: x.imag, y1_3)), where='post', label = 'Im')
ax1_3.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y1_3)), where='post', label = 'Mod')
ax1_3.set_xlabel("Tn (repetition period number)")
ax1_3.set_ylabel("Signal")
ax1_3.legend()
ax1_3.grid(True)
# Спектр азимутальной пачки
fig1_3_2, ax1_3_2 = plt.subplots()
ax1_3_2.set_title('Spectrum of Azimuth graph (LFM signal), target 3')
ax1_3_2.step(x1, abs(np.fft.fft(y1_3)))
ax1_3_2.legend()
ax1_3_2.grid(True)
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
# Первая цель
fig3_1, ax3_1 = plt.subplots()
ax3_1.set_title('Azimuth graph (simple rectangular modulation), target 1')
ax3_1.step(x1, list(map(lambda x: x.real, y3_1)), where='post', label = 'Re')
ax3_1.step(x1, list(map(lambda x: x.imag, y3_1)), where='post', label = 'Im')
ax3_1.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y3_1)), where='post', label = 'Mod')
ax3_1.set_xlabel("Tn (repetition period number)")
ax3_1.set_ylabel("Signal")
ax3_1.legend()
ax3_1.grid(True)
# Спектр азимутальной пачки
fig3_1_2, ax3_1_2 = plt.subplots()
ax3_1_2.set_title('Spectrum of Azimuth graph (simple rectangular modulation), target 1')
ax3_1_2.step(x1, abs(np.fft.fft(y3_1)))
ax3_1_2.legend()
ax3_1_2.grid(True)
# Вторая цель
fig3_2, ax3_2 = plt.subplots()
ax3_2.set_title('Azimuth graph (simple rectangular modulation), target 2')
ax3_2.step(x1, list(map(lambda x: x.real, y3_2)), where='post', label = 'Re')
ax3_2.step(x1, list(map(lambda x: x.imag, y3_2)), where='post', label = 'Im')
ax3_2.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y3_2)), where='post', label = 'Mod')
ax3_2.set_xlabel("Tn (repetition period number)")
ax3_2.set_ylabel("Signal")
ax3_2.legend()
ax3_2.grid(True)
# Спектр азимутальной пачки
fig3_2_2, ax3_2_2 = plt.subplots()
ax3_2_2.set_title('Spectrum of Azimuth graph (simple rectangular modulation), target 2')
ax3_2_2.step(x1, abs(np.fft.fft(y3_2)))
ax3_2_2.legend()
ax3_2_2.grid(True)
# Третья цель
fig3_3, ax3_3 = plt.subplots()
ax3_3.set_title('Azimuth graph (simple rectangular modulation), target 3')
ax3_3.step(x1, list(map(lambda x: x.real, y3_3)), where='post', label = 'Re')
ax3_3.step(x1, list(map(lambda x: x.imag, y3_3)), where='post', label = 'Im')
ax3_3.step(x1, list(map(lambda x: np.sqrt(np.real(x)**2 + np.imag(x)**2), y3_3)), where='post', label = 'Mod')
ax3_3.set_xlabel("Tn (repetition period number)")
ax3_3.set_ylabel("Signal")
ax3_3.legend()
ax3_3.grid(True)
# Спектр азимутальной пачки
fig3_3_2, ax3_3_2 = plt.subplots()
ax3_3_2.set_title('Spectrum of Azimuth graph (simple rectangular modulation), target 3')
ax3_3_2.step(x1, abs(np.fft.fft(y3_3)))
ax3_3_2.legend()
ax3_3_2.grid(True)
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