import numpy as np
import csv

# Задание параметров радара
class Radar:

	def __init__(self, Tn, ti, Tws, fs, Fz, modLaw, disp, speed):
		self.Tn = Tn
		self.ti = ti
		self.Tws = Tws
		self.fs = fs
		self.Fz = Fz
		self.disp = 1
		self.dR = int(300000000 / (2 * fs))
		self.Rmax = int(300000000 * ( Tws + ti)/self.dR/2)
		print ("Rmax = " + str(self.Rmax))
		self.modLaw = np.zeros(self.Rmax, dtype = np.complex)
		self.speed = speed
		self.bettaA = Tn*speed
		if(modLaw == "Rect"):
			for i in range(len(modLaw)):
				self.modLaw[i] = 1 if i <= ti else 0

# Задание параметров цели
class Target:

	def __init__(self, snr, Vr, D, R, tr, B):
		self.snr = snr
		self.Vr = Vr
		self.D = D
		self.Dfin = D
		self.R = R
		self.tr = tr
		self.B = B
		self.Bfin = B

# Задание параметров мешающего отражения
class SClutter(Target):

	def __init__(self, snr, Vr, D_start, D_finish, R, tr, B_start, B_finish):
		super().__init__(snr, Vr, D_start, R, tr, B_start)
		self.Dfin = D_finish
		self.Bfin = B_finish

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
			Amp = 10**(self.target.snr/20)
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

# Имитация отражённого от целей сигнала используя зарание заданные параметры целей и РЛС
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
		for i in range(count):
			Re = np.random.normal(0, self.radar.disp)
			Im = np.random.normal(0, self.radar.disp)
			yield np.complex(Re,Im)

	# Сохранение результатов расчёта развёртки дальности для каждого периода повторений в файл
	def saveRes(sweepRange):
		with open("sweepRange.csv", "a") as f:
			writer = csv.writer(f)
			writer.writerow(sweepRange)