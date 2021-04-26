import numpy as np
import csv


"""
Класс, задающий параметры радара.

Входные аргументы констуктора:
	Tn - периол повторения
	ti - длительность импульса
	Tws - период приёма отражённого сигнала
	fs - частота дискретизации
	Fz - частота зондирующего сигнала
	P - мощность передатчика
	Gt - коэффициент усиления передатчика
	Gr - коэффициент усиления приемника
	dist - дисперсия внутреннего шума приёмника
	dR - дискрет дальности
	Rmax - максимальная дальность работы
	modLaw - закон модуляции
	speed - скорость вращения антенны
	bettaA - начальный угол азимута антены
"""
class Radar:

	def __init__(self, Tn, ti, Tws, fs, Fz, modLaw, disp, speed, P, Gt):
		self.Tn = Tn
		self.ti = ti
		self.Tws = Tws
		self.fs = fs
		self.Fz = Fz
		self.P = P
		self.Gt = Gt
		self.Gr = 40000/70
		self.disp = 1
		self.dR = int(300000000 / (2 * fs))
		self.Rmax = int(300000000 * ( Tws + ti)/self.dR/2)
		self.modLaw = np.zeros(self.Rmax, dtype = np.complex)
		self.speed = speed
		self.bettaA = Tn*speed
		if(modLaw == "Rect"):
			for i in range(len(modLaw)):
				self.modLaw[i] = np.complex(1, 0) if i <= ti else np.complex(0, 0)
		if(modLaw == "LFM"):
			Fdev = 750000
			self.Fdev = Fdev
			length = int(ti*fs+0.5)
			dt = 1/fs
			for i in range(len(modLaw)):
				if i in range(length):
					t = -ti/2+i*dt
					frequency = 2 * np.pi * (abs(Fdev)*(t-t*t)/(2*ti))
					self.modLaw[i] = np.complex(np.cos(frequency), np.sin(frequency))
				else:
					self.modLaw[i] = np.complex(0, 0)

"""
Класс, задающий параметры цели.

Входные аргументы констуктора:
	snr - соотношение сигнал/шум
	Vr - радиальная скорость
	D - дальность
	R - закон корреляции
	tr - время корреляции
	B - угол азимута
	esa - эффективная площадь рассеивания
"""
class Target:

	def __init__(self, snr, Vr, D, R, tr, B, esa):
		self.snr = snr
		self.esa = esa
		self.Vr = Vr
		self.D = D
		self.Dfin = D
		self.R = R
		self.tr = tr
		self.B = B
		self.Bfin = B

"""
Класс, задающий параметры для мешающего отражения

Входные аргументы констуктора:
	snr - соотношение сигнал/шум
	Vr - радиальная скорость
	D_start - дальность до начала
	D_finish - дальность до конца
	R - закон корреляции
	tr - время корреляции
	B_start - угол азимута начала
	B_finish - угол азимута конца
	esa - эффективная площадь рассеивания
"""
class SClutter(Target):

	def __init__(self, snr, Vr, D_start, D_finish, R, tr, B_start, B_finish, esa):
		super().__init__(snr, Vr, D_start, R, tr, B_start, esa)
		self.Dfin = D_finish
		self.Bfin = B_finish

"""
Класс, выполняющий рассчёт отражённого сигнала от цели для каждого периода повторений.

Входные аргументы констуктора:
	target - объект класса Target
	radar - объект класса Radar
	fi1ter - объект класса Filter
	count - число периодов повторения
"""
class CalcTargetSignal():

	def __init__(self, target, radar, fi1ter, count):
		self.target = target
		self.radar = radar
		self.fi1ter = fi1ter
		self.dR = radar.dR
		self.targetfi1ter = list()
		self.maxEl = 0
		self.targetR = self.targetFilter(count)

	# Функция возвращающая генератор отражённого сигнала от цели
	def targetFilter(self, count):
		def targetGen(i=0):
			size = int((self.target.Dfin - self.target.D)/self.dR) + 1
			buffer = list(map(lambda x: x*2/self.maxEl , self.targetfi1ter[-(1+i):-(size+i+1):-1]))
			i += size
			yield buffer
		if(self.target.R == "gaussian"):
			self.targetfi1ter = list(self.fi1ter.gaussian(1000 + count*(int((self.target.Dfin - self.target.D)/self.dR) + 1)))
			self.maxEl = max(self.targetfi1ter)
			return targetGen
		if(self.target.R == "exponential"):
			self.targetfi1ter = list(self.fi1ter.exponential(1000 + count*(int((self.target.Dfin - self.target.D)/self.dR) + 1)))
			self.maxEl = max(self.targetfi1ter)
			return targetGen
		if(self.target.R == "exponent_parabolic_1"):
			self.targetfi1ter = list(self.fi1ter.exponent_parabolic_1(1000 + count*(int((self.target.Dfin - self.target.D)/self.dR) + 1)))
			self.maxEl = max(self.targetfi1ter)
			return targetGen
		if(self.target.R == "exponent_parabolic_2"):
			self.targetfi1ter = list(self.fi1ter.exponent_parabolic_2(1000 + count*(int((self.target.Dfin - self.target.D)/self.dR) + 1)))
			self.maxEl = max(self.targetfi1ter)
			return targetGen

	# Функция рассчитывающая отражённый сигнал от цели учитывая поправку Доплера и ДНА
	def calc(self, numTn, isClut = False):
		listFilter = list()
		targetSignalBuff = next(self.targetR())
		doplerPhase = self.doplerPhase(numTn)
		for i in range(len(targetSignalBuff)):
			a = np.real(targetSignalBuff[i])
			b = np.imag(targetSignalBuff[i])
			mod = np.sqrt(a**2 + b**2)
			F = 0
			if a > 0 and b > 0 : F = 0
			if a < 0 and b > 0 : F = np.pi
			if a < 0 and b < 0 : F = -np.pi
			if a > 0 and b < 0 : F = 0
			currentPhase = np.arctan(b/a) + F
			totalPhase = currentPhase + doplerPhase
			Re = mod*np.cos(totalPhase)
			Im = mod*np.sin(totalPhase)
			# Расчёт амплитуды через соотношения сигнал/шум:
			#Amp = 10**(self.target.snr/20)
			P = self.radar.P*self.radar.Gt*self.radar.Gr*self.target.esa*(300000000/self.radar.Fz)/((4*np.pi)**3*(self.target.D+i)**4)
			Amp = np.sqrt(P)
			self.lastEl = np.complex(Re, Im)
			listFilter.append(Amp * self.getAP(numTn) * targetSignalBuff[i])
		return listFilter

	# Функция расчёта поправки Доплера для каждого периода повторений
	def doplerPhase(self, numTn):
		l = 300000000 / self.radar.Fz
		doplerPhase = 4*np.pi*self.target.Vr/l*self.radar.Tn*numTn
		return doplerPhase

	# Функция возвращающая ДНА для каждого периода повторений
	def getAP(self, numTn):
		bettaA = self.radar.bettaA*numTn % 360
		Bstart = self.target.B
		Bfin = self.target.Bfin
		if Bstart > Bfin:
			if (bettaA > Bstart and bettaA > Bfin) or (bettaA < Bstart and bettaA < Bfin): return 1
		if ((bettaA < Bstart) or (bettaA > Bfin)) and (bettaA != Bstart and bettaA != Bfin):
			Amin = min(abs(bettaA - Bstart), abs(bettaA - Bfin))
			return abs(np.sin(Amin)/Amin)
		else:
			return 1

"""
Класс, имитирующий отражённый от целей сигнал, используя зарание заданные параметры целей и РЛС.

Входные аргументы констуктора:
	radar - объект класса Radar
	targets - список объектов класса Target
	fi1ter - объект класса Filter
"""
class ImSignal:

	def __init__(self, radar, targets, fi1ter):
		self.radar = radar
		self.targets = targets
		self.fi1ter = fi1ter
		self.sweepRange = np.zeros(self.radar.Rmax, dtype = np.complex)

	# Основная функция расчёта развёртки дальности, для каждого периода повторений.
	def main(self, count):
		setSweepRange = list()
		calcTargetSignal = list()
		for target in self.targets:
			calcTargetSignal.append(CalcTargetSignal(target, self.radar, self.fi1ter, count))
		for k in range(1, count+1):
			self.sweepRange = np.zeros(self.radar.Rmax, dtype = np.complex)
			for obj in calcTargetSignal:
				listCalc = obj.calc(k)
				for i in range(len(listCalc)):
					self.sweepRange[int(i + obj.target.D/self.radar.dR)] += listCalc[i]
			ImSignal.compress(self.sweepRange, self.radar.modLaw)
			for j, i in enumerate(self.coldNoiseGen(self.radar.Rmax)):
				self.sweepRange[j] += i
			ImSignal.saveRes(self.sweepRange)
			setSweepRange.append(self.sweepRange)
		return setSweepRange

	# Функция сжатия развёртки дальности и закона модуляции РЛС
	def compress(sweepRange, modLaw):
		fftSweepRange = np.fft.fft(sweepRange)
		fftModLaw = np.fft.fft(modLaw)
		return np.fft.ifft(fftSweepRange*fftModLaw)

	# Генератор внутренних шумов принимающего устройства РЛС
	def coldNoiseGen(self, count):
		k = 1.38*10**(-23)
		T = 300
		Bn = 4
		Amp = Bn*T*k
		for i in range(count):
			# при расчёте через соотношение сигнал/шум
			#------------------------------------------
			#Re = np.random.normal(0, self.radar.disp)
			#Im = np.random.normal(0, self.radar.disp)
			#------------------------------------------
			# при расчёте через основную формулу радиолокации
			#------------------------------------------
			Re = np.random.normal(0, Amp)
			Im = np.random.normal(0, Amp)
			#------------------------------------------
			yield np.complex(Re,Im)

	# Сохранение результатов расчёта развёртки дальности для каждого периода повторений в файл
	def saveRes(sweepRange):
		with open("sweepRange.csv", "a") as f:
			writer = csv.writer(f)
			writer.writerow(sweepRange)