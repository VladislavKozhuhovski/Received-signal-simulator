from matplotlib import pyplot as plt
from myfilter import Filter
from fourier import Fourier_transform

# Задание фильтров
filter1 = Filter(0.003, 0.1, isNormal = True)
filter2 = Filter(0.003, 0.1, isNormal = False)

# Реализация фильтров для нормального распределения
yGauss = list(filter1.gaussian(2048))
yExp = list(filter1.exponential(2048))
yExpParab1 = list(filter1.exponent_parabolic_1(2048))
yExpParab2 = list(filter1.exponent_parabolic_2(2048))

# Построение графиков реализации фильтров
fig1, ax1 = plt.subplots()
ax1.set_title('Gaussian and exponential filter for normal distribution')
ax1.plot(yGauss, label="Gaussian filter")
ax1.plot(yExp, label="Exponential filter")
ax1.legend()
ax1.grid(True)
fig2, ax2 = plt.subplots()
ax2.set_title('Exponent-parabolic filter for normal distribution')
ax2.plot(yExpParab1)
ax2.grid(True)
fig3, ax3 = plt.subplots()
ax3.set_title('Exponent-parabolic filter for normal distribution')
ax3.plot(yExpParab2)
ax3.grid(True)

# Спектр для нормального распределения
YGauss = abs(Fourier_transform().FFT(yGauss))
YExp = abs(Fourier_transform().FFT(yExp))
YExpParab1 = abs(Fourier_transform().FFT(yExpParab1))
YExpParab2 = abs(Fourier_transform().FFT(yExpParab2))

# Построение графиков спектра для нормального распределения
xAxis = [i for i in range(1048)]
fig4, ax4 = plt.subplots()
ax4.set_title('Spectrum of gaussian and exponential filter for normal distribution')
ax4.bar(xAxis, YGauss[1000:], label="Gaussian filter")
ax4.bar(xAxis, YExp[1000:], label="Exponential filter")
ax4.legend()
ax4.grid(True)
fig5, ax5 = plt.subplots()
ax5.set_title('Spectrum of Exponent-parabolic filter for normal distribution')
ax5.bar(xAxis, YExpParab1[1000:])
ax5.grid(True)
fig6, ax6 = plt.subplots()
ax6.set_title('Spectrum of exponent-parabolic filter for normal distribution')
ax6.bar(xAxis, YExpParab2[1000:])
ax6.grid(True)

# Реализация фильтров для равномерного распределения
zGauss = list(filter2.gaussian(2048))
zExp = list(filter2.exponential(2048))
zExpParab1 = list(filter2.exponent_parabolic_1(2048))
zExpParab2 = list(filter2.exponent_parabolic_2(2048))

# Построение графиков реализации фильтров
fig7, ax7 = plt.subplots()
ax7.set_title('Gaussian and exponential filter for uniform distribution')
ax7.plot(zGauss, label="Gaussian filter")
ax7.plot(zExp, label="Exponential filter")
ax7.legend()
ax7.grid(True)
fig8, ax8 = plt.subplots()
ax8.set_title('Exponent-parabolic filter for uniform distribution')
ax8.plot(zExpParab1)
ax8.grid(True)
fig9, ax9 = plt.subplots()
ax9.set_title('Exponent-parabolic filter for uniform distribution')
ax9.plot(zExpParab2)
ax9.grid(True)

# Спектр для равномерного распределения
ZGauss = abs(Fourier_transform().FFT(zGauss))
ZExp = abs(Fourier_transform().FFT(zExp))
ZExpParab1 = abs(Fourier_transform().FFT(zExpParab1))
ZExpParab2 = abs(Fourier_transform().FFT(zExpParab2))

# Построение графиков спектра для нормального распределения
xAxis = [i for i in range(1048)]
fig10, ax10 = plt.subplots()
ax10.set_title('Spectrum of gaussian and exponential filter for uniform distribution')
ax10.bar(xAxis, ZGauss[1000:], label="Gaussian filter")
ax10.bar(xAxis, ZExp[1000:], label="Exponential filter")
ax10.legend()
ax10.grid(True)
fig11, ax11 = plt.subplots()
ax11.set_title('Spectrum of exponent-parabolic filter for uniform distribution')
ax11.bar(xAxis, ZExpParab1[1000:])
ax11.grid(True)
fig12, ax12 = plt.subplots()
ax12.set_title('Spectrum of exponent-parabolic filter for uniform distribution')
ax12.bar(xAxis, ZExpParab2[1000:])
ax12.grid(True)

plt.show()