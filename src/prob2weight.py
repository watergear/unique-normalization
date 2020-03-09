from enum import Enum
import math

W = 256
WE = 20
ScatterType = Enum('ScatterType', ('Smooth', 'AccurateBoundary'))
scatter_type = ScatterType.Smooth

class Prob2Weight:
	def __init__(self, K, A, p0_K):
		self.W = W
		self.WE = WE
		self.scatter_type = scatter_type

		self.K = K
		self.A = A
		self.p0 = p0_K / K

		self.L = (self.W - 2*self.WE) / (2*math.log(self.A))

		delta = 0.5
		self.D = -1 * math.log(math.exp(delta/self.L) - 1)
		if self.scatter_type is ScatterType.Smooth:
			self.tk = ((self.W-2)/(2*self.L) - math.log(self.p0)) / math.log(self.K) # smooth
		else:
			self.tk = (math.log(self.A) - math.log(self.p0) + self.D) / math.log(self.K) # accurate boundary

		self._K_tk = math.pow(self.K, self.tk)

		self._h_p0 = self.h(self.p0)

		self._p0_A_ = self.p0/self.A
		self._p0_A = self.p0*self.A
		self.G0 = self.g(self._p0_A_)
		self.G1 = self.g(self._p0_A)

		self._W_1_G0_g_0_G0 = (self.W-1-self.G0)/(self.g(0)-self.G0)
		self._1_G1_g_1_G1 = (1-self.G1)/(self.g(1)-self.G1)

	def Weight(self, p):
		return int(self.f(p)+0.5)

	# Weight(x*m) = Weight(x) + MultipleWeight(m)
	def MultipleWeight(self, m):
		return int(self.f_m(m)+0.5)

	def h(self, x):
		return math.log(x*self._K_tk + 1)
	def g(self, x):
		return (self._h_p0 - self.h(x))*self.L + 0.5*self.W
	def f(self, x):
		g_x = self.g(x)
		if x < self._p0_A_:
			g_x = self.G0 + (g_x-self.G0)*self._W_1_G0_g_0_G0
		elif self._p0_A < x:
			g_x = self.G1 + (g_x-self.G1)*self._1_G1_g_1_G1
		return g_x
	def f_m(self, x):
		return -1*math.log(x)*self.L

def test():
	K = 141235
	a = 4
	ta = 4

	p2w = Prob2Weight(K, math.pow(a, ta), 1)

	print("K = %d" %(p2w.K))
	print("W = %d" %(p2w.W))
	print("WE = %d" %(p2w.WE))
	print("A = %f" %(p2w.A))
	print("p0 = %e" %(p2w.p0))
	print("L = %f" %(p2w.L))
	print("G0 = %f" %(p2w.G0))
	print("G1 = %f" %(p2w.G1))
	print("scatter type = %s" %("smooth" if ScatterType.Smooth == p2w.scatter_type else "accurate boundary"))
	# print("scatter type = %s" %(p2w.scatter_type))
	print("D = %f" %(p2w.D))
	print("tk = %f" %(p2w.tk))

	print("")
	print("prob %d -> origin weight %f" %(0, p2w.g(0)))
	print("prob %d -> origin weight %f" %(1, p2w.g(1)))
	print("prob %d -> weight %f" %(0, p2w.f(0)))
	print("prob %d -> weight %f" %(1, p2w.f(1)))

	print("")
	probs = []
	for i in range(-1*ta, ta+1):
		probs.append(math.pow(a, i) * p2w.p0)
	for p in probs:
		print("prob %e -> weight %f" %(p, p2w.f(p)))
	print("")
	prob_multis = []
	for i in range(1, 10+1):
		prob_multis.append(i*0.1)
	for m in prob_multis:
		print("prob multi %f -> weight %f" %(m, p2w.f_m(m)))

def demo():
	K = 141235
	a = 4
	ta = 4
	p2w = Prob2Weight(K, math.pow(a, ta), 1)
	probs = [0, 0.01/K, 0.1/K, 1.0/K, 10.0/K, 100.0/K, 1.0/2, 1]
	for p in probs:
		print("prob %e -> weight %d" %(p, p2w.Weight(p)))
	prob_multis = [0.25, 0.5, 0.75, 1]
	print("")
	for m in prob_multis:
		print("prob multi %f -> weight %d" %(m, p2w.MultipleWeight(m)))

if __name__ == "__main__":
	# execute only if run as a script
	print("-------- test --------")
	test()

	print("-------- demo --------")
	demo()