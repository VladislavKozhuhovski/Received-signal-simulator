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
		self.R = R
		self.tr = tr
		self.B = B

class SClutter(Target):

	def __init__(self, snr, Vr, D_start, D_finish, R, tr, B_start, B_finish):
		super().__init__(snr, Vr, D_start, R, tr, B_start)
		self.Dfin = D_finish
		self.Bfin = B_finish


class CalcTargetSignal():

	def __init__(self, target, radar, fi1ter, numTn):
		self.target = target
		self.radar = radar
		self.fi1ter = fi1ter
		self.numTn = numTn

	def calc(self):
		if(self.target.R == "gaussian"):
			gausfi1ter = list(self.fi1ter.gaussian(1000))
			maxEl = max(gausfi1ter)
			lastEl = gausfi1ter[-1]/(2*maxEl)
		if(self.target.R == "exponential"):
			expfi1ter = list(self.fi1ter.exponential(1000))
			maxEl = max(expfi1ter)
			lastEl = expfi1ter[-1]/(2*maxEl)
		if(self.target.R == "exponent_parabolic_1"):
			exp_parab_fi1ter = list(self.fi1ter.exponent_parabolic_1(1000))
			maxEl = max(exp_parab_fi1ter)
			lastEl = exp_parab_fi1ter[-1]/(2*maxEl)
		if(self.target.R == "exponent_parabolic_2"):
			exp_parab_fi1ter = list(self.fi1ter.exponent_parabolic_2(1000))
			maxEl = max(exp_parab_fi1ter)
			lastEl = exp_parab_fi1ter[-1]/(2*maxEl)
		a = np.real(lastEl)
		b = np.imag(lastEl)
		mod = np.sqrt(a**2 + b**2)
		F = 0
		if a > 0 and b > 0 : F = 0
		if a < 0 and b > 0 : F = np.pi
		if a < 0 and b < 0 : F = -np.pi
		if a > 0 and b < 0 : F = 0
		currentPhase = np.arctan(b/a) + F
		l = 300000000 / self.radar.Fz
		doplerPhase = 4*np.pi*self.target.Vr/l*self.radar.Tn*self.numTn
		totalPhase = currentPhase + doplerPhase
		Re = mod*np.cos(totalPhase)
		Im = mod*np.sin(totalPhase)
		Amp = 10**(self.target.snr/20)
		lastEl = np.complex(Re, Im)
		return Amp * self.getAP() * lastEl

	def getAP(self):
		bettaA = self.radar.bettaA*self.numTn
		bettaT = self.target.B
		return abs(np.sin(bettaA - bettaT)/(bettaA - bettaT))

# Имитация отражённого от целей сигнала используя зарание заданные параметры целей и РЛС
class ImSignal:

	def __init__(self, radar, targets, fi1ter):
		self.radar = radar
		self.targets = targets
		self.fi1ter = fi1ter
		self.sweepRange = np.zeros(4096, dtype = np.complex)


	def main(self, k):
		self.sweepRange = np.zeros(4096, dtype = np.complex)
		dR = int(300000000 / (2 * self.radar.fs))
		for target in self.targets:
			if hasattr(target, 'Dfin'):
				for i in range(0, int((target.Dfin - target.D)/dR)):
					self.sweepRange[int((i+target.D)/dR)] = CalcTargetSignal(target, self.radar, self.fi1ter, k).calc()
			else:
				self.sweepRange[int(target.D/dR)] = CalcTargetSignal(target, self.radar, self.fi1ter, k).calc()
		ImSignal.compress(self.sweepRange, self.radar.modLaw)
		for j, i in enumerate(self.coldNoiseGen(4096)):
			self.sweepRange[j] += i
		ImSignal.saveRes(self.sweepRange)
		return self.sweepRange

	def compress(sweepRange, modLaw):
		fftSweepRange = np.fft.fft(sweepRange)
		fftModLaw = np.fft.fft(modLaw)
		print(fftSweepRange*fftModLaw)
		return np.fft.rfft(fftSweepRange*fftModLaw)

	def coldNoiseGen(self, count):
		for i in range(count):
			Re = np.random.normal(0, self.radar.disp)
			Im = np.random.normal(0, self.radar.disp)
			yield np.complex(Re,Im)

	def saveRes(sweepRange):
		with open("sweepRange.csv", "a") as f:
			writer = csv.writer(f)
			writer.writerow(sweepRange)