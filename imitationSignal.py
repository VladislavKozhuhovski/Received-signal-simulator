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
		#self.Tc = Tc
		self.modLaw = np.zeros(4096, dtype = np.complex)
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

class SClutter(Target):

	def __init__(self, snr, Vr, D_start, D_finish, R, tr, B_start, B_finish):
		super().__init__(snr, Vr, D_start, R, tr, B_start)
		self.Dfin = D_finish
		self.Bfin = B_finish


class CalcTargetSignal():

	def __init__(self, target, radar, fi1ter):
		self.target = target
		self.radar = radar
		self.fi1ter = fi1ter
		self.lastEl = self.targetFilter()

	def targetFilter(self):
		if(self.target.R == "gaussian"):
			gausfi1ter = list(self.fi1ter.gaussian(1000))
			maxEl = max(gausfi1ter)
			lastEl = gausfi1ter[-1]*2/(maxEl)
		if(self.target.R == "exponential"):
			expfi1ter = list(self.fi1ter.exponential(1000))
			maxEl = max(expfi1ter)
			lastEl = expfi1ter[-1]*2/(maxEl)
		if(self.target.R == "exponent_parabolic_1"):
			exp_parab_fi1ter = list(self.fi1ter.exponent_parabolic_1(1000))
			maxEl = max(exp_parab_fi1ter)
			lastEl = exp_parab_fi1ter[-1]*2/(maxEl)
		if(self.target.R == "exponent_parabolic_2"):
			exp_parab_fi1ter = list(self.fi1ter.exponent_parabolic_2(1000))
			maxEl = max(exp_parab_fi1ter)
			lastEl = exp_parab_fi1ter[-1]*2/(maxEl)
		return lastEl

	def calc(self, numTn, isClut = False):
		listFilter = list()
		dR = int(300000000 / (2 * self.radar.fs))
		for i in range(int((self.target.Dfin - self.target.D)/dR) + 1):
			a = np.real(self.lastEl)
			b = np.imag(self.lastEl)
			mod = np.sqrt(a**2 + b**2)
			F = 0
			if a > 0 and b > 0 : F = 0
			if a < 0 and b > 0 : F = np.pi
			if a < 0 and b < 0 : F = -np.pi
			if a > 0 and b < 0 : F = 0
			currentPhase = np.arctan(b/a) + F
			l = 300000000 / self.radar.Fz
			doplerPhase = 4*np.pi*self.target.Vr/l*self.radar.Tn*numTn
			totalPhase = currentPhase + doplerPhase
			Re = mod*np.cos(totalPhase)
			Im = mod*np.sin(totalPhase)
			Amp = 10**(self.target.snr/20)
			self.lastEl = np.complex(Re, Im)
			listFilter.append(Amp * self.getAP(numTn) * self.lastEl)
		return listFilter

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
		self.sweepRange = np.zeros(4096, dtype = np.complex)


	def main(self, count):
		dR = int(300000000 / (2 * self.radar.fs))
		setSweepRange = list()
		for k in range(1, count+1):
			self.sweepRange = np.zeros(4096, dtype = np.complex)
			calcTargetSignal = list()
			for target in self.targets:
				calcTargetSignal.append(CalcTargetSignal(target, self.radar, self.fi1ter))
			for obj in calcTargetSignal:
				listCalc = obj.calc(k)
				for i in range(len(listCalc)):
					self.sweepRange[int(i + obj.target.D/dR)] += listCalc[i]
			ImSignal.compress(self.sweepRange, self.radar.modLaw)
			for j, i in enumerate(self.coldNoiseGen(4096)):
				self.sweepRange[j] += i
			ImSignal.saveRes(self.sweepRange)
			setSweepRange.append(self.sweepRange)
		return setSweepRange


	def compress(sweepRange, modLaw):
		fftSweepRange = np.fft.fft(sweepRange)
		fftModLaw = np.fft.fft(modLaw)
		return np.fft.ifft(fftSweepRange*fftModLaw)

	def coldNoiseGen(self, count):
		for i in range(count):
			Re = np.random.normal(0, self.radar.disp)
			Im = np.random.normal(0, self.radar.disp)
			yield np.complex(Re,Im)

	def saveRes(sweepRange):
		with open("sweepRange.csv", "a") as f:
			writer = csv.writer(f)
			writer.writerow(sweepRange)