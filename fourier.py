import numpy as np
class Fourier_transform:

	def __init__(self):
		pass

	# ПФ
	def DFT_slow(self, x):
		"""Compute the discrete Fourier Transform of the 1D array x"""
		x = np.asarray(x, dtype=float)
		N = x.shape[0]
		n = np.arange(N)
		k = n.reshape((N, 1))
		M = np.exp(-2j * np.pi * k * n / N)
		return np.dot(M, x)

	# БПФ
	def FFT(self, x):
		"""A recursive implementation of the 1D Cooley-Tukey FFT"""
		x = np.asarray(x, dtype=float)
		N = x.shape[0]
		#print("N = " + str(N))
		if N % 2 > 0:
			raise ValueError("size of x must be a power of 2")
		elif N <= 32:
			return self.DFT_slow(x)
		else:
			X_even = self.FFT(x[::2])
			X_odd = self.FFT(x[1::2])
			factor = np.exp(-2j * np.pi * np.arange(N) / N)
			return np.concatenate([X_even + factor[: int(N / 2)] * X_odd, X_even + factor[int(N / 2):] * X_odd])