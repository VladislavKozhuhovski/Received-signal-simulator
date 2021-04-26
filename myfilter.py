import numpy as np
class Filter:

	# N - объём выборки сигнала
	def __init__(self, Tn, t, isNormal = False):
		self.Tn = Tn
		self.isNormal = isNormal
		self.t = t
		self.mean = 1
		self.disp = 1

	# генерация случайного значения по нормальнаму закону распределения
	def _normal_distribution(self):
		Re = np.random.normal(0, 1/self.Tn)
		Im = np.random.normal(0, 1/self.Tn)
		x = np.complex(Re,Im)
		return x

	# генерация случайного значения по равномерному закону распределения
	def _uniform_distribution(self):
		mod = np.random.uniform(0, 1)
		fi = np.random.uniform(0, 2*np.pi)
		x = np.complex(mod*np.cos(fi), mod*np.sin(fi))
		return x

	# Фильтр Гаусса
	def gaussian(self, count):
		if count == None: count = 1
		n = 500
		sumCk, p = 0, 0
		alfa = 2/self.t
		error = 0.001
		gamma = alfa*self.Tn
		if(gamma > 0.5):
			raise ValueError("gamma for gaussian filter is not correct")
		for j in range(n):
			ck = np.zeros(2*j+1)
			sumCk = 0
			for i in range(2*j+1):
				ck[i] = 2*self.mean*(np.sqrt(gamma/np.pi))*1/(1+4*(gamma**2)*((i-j)**2))
				#ck[i] = self.mean*np.sqrt(2*gamma)/np.pi**(1/4)*np.exp(-2*gamma**2*i**2)
				#ck[i] = self.mean/(np.sqrt(np.pi*gamma)*np.sin(gamma*i)/i)
				sumCk += ck[i]**2
			if (abs(1-1/self.disp*sumCk) >= error): 
				#print("Error = " + str(abs(1-1/self.disp*sumCk)))
				if j == n-1: raise ValueError("Gaussian filter couldn't find the correct p")
				continue
			else: 
				p = j
				#print("p = " + str(j))
				break
		if self.isNormal == True: x = [self._normal_distribution() for i in range(count + p)]
		else: x = [self._uniform_distribution() for i in range(count + p)]
		for i in range(count):
			y = 0
			for j in range(-p,p+1):
				if i-j >= 0 and i-j < count:
					y += ck[j+p]*x[i-j]
				else:
					continue
			yield y

	# Экспоненциальный фильтр
	def exponential(self, count):
		if count == None: count = 1
		y = np.zeros(count, dtype=np.complex)
		a0 = np.sqrt(2*self.t)*self.Tn/(self.Tn+self.t)
		b0 = self.t/(self.Tn+self.t)
		if self.isNormal == True: x = [self._normal_distribution() for i in range(count)]
		else: x = [self._uniform_distribution() for i in range(count)]
		for i in range(count):
			if(i==0):
				y[i] = a0*x[i]
			else:
				y[i] = a0*x[i] + b0*y[i-1]
			yield y[i]

	# Экспоненциально-параболический фильтр первая реализация
	def exponent_parabolic_1(self, count):
		if count == None: count = 1
		res = np.zeros(count, dtype=np.complex)
		alfa = 2/self.t
		gamma = alfa*self.Tn
		p = np.exp(-gamma)
		alfa1 = 1-4*(p**2)*gamma-p**4
		alfa0 = (p**3)*(1+gamma)-p*(1+gamma)
		b2 = -(p**2)
		b1 = 2*p
		a0 = self.mean*np.sqrt((alfa1**2 + np.sqrt((alfa1**2) - 4*(alfa0**2), dtype=np.complex))/2)
		a1 = self.mean*alfa0/alfa1
		if self.isNormal == True: x = [self._normal_distribution() for i in range(count)]
		else: x = [self._uniform_distribution() for i in range(count)]
		res[0] = a0*x[0]
		yield res[0]
		res[1] = a0*x[1] + a1*x[0] + b1*res[0]
		yield res[1]
		for i in range(2, count):
			res[i] = a0*x[i] + a1*x[i-1]+b1*res[i-1]+b2*res[i-2]
			yield res[i]

	# Экспоненциально-параболический фильтр вторая реализация
	def exponent_parabolic_2(self, count):
		if count == None: count = 1
		res = np.zeros(count, dtype=np.complex)
		p = 2*self.Tn+self.t
		b2 = -self.t**2/p**2
		b1 = 2*self.t/p
		a0 = 4*np.sqrt(2*self.t)*self.Tn**2/(p**2)
		if self.isNormal == True: x = [self._normal_distribution() for i in range(count)]
		else: x = [self._uniform_distribution() for i in range(count)]
		res[0] = a0*x[0]
		yield res[0]
		res[1] = a0*x[1] + b1*res[0]
		yield res[1]
		for i in range(2, count):
			res[i] = a0*x[i]+b1*res[i-1]+b2*res[i-2]
			yield res[i]