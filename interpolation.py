import matplotlib.pyplot as plt
import numpy as np
import time

yp=np.array([139620, 173706, 311622])
xp=np.array([1950, 1960, 1970])

ypo=np.array([311622, 633392, 1011500])
xpo=np.array([1970, 1980, 1990])

ypq=np.array([139620, 173706, 311622, 633392])
xpq=np.array([1950, 1960, 1970, 1980])

ypoq=np.array([311622, 633392, 1011500, 1403793])
xpoq=np.array([1970, 1980, 1990, 2000])

xnp = np.array([0.0, 0.2, 0.4, 0.6, 0.8])
ynp = np.array([1.008, 1.064, 1.125, 1.343, 1.512])

#Parent Method
class InterpolationMethod():
    Xint, Yint = None, None
    time = 0

    def __init__(self, x, y):
        self.setInput(x, y)

    def setInput(self, x, y):
        self.Xint = x
        self.Yint = y

    def graph(self, function, legend, xlimit, ylimit, dots):
        t = np.linspace(xlimit, ylimit, dots)
        y = np.zeros(len(t))

        for i in range(len(t)):
            y[i] = function(t[i])
        plt.plot(t, y, label=legend)
        plt.plot(self.Xint, self.Yint, 'bo')
        plt.plot(self.Xint, self.Yint, 'bo', label = 'Pontos inseridos')
        plt.title(legend)
        plt.legend();
        plt.ylim([0, self.Yint[self.Yint.shape[0]-1]*1.4])
        plt.show()

class Polinear(InterpolationMethod):
    def __init__(self, x, y):
        InterpolationMethod.__init__(self,x,y)

    def generatePolynomial(self, M):
        matx = np.zeros((M.shape[0], M.shape[0]))
        for i in range (0, M.shape[0]):
            for k in range(0, M.shape[0]):
                matx[i][k] = pow(M[i],M.shape[0]-(k+1))
        return matx

    def f(self, t, h, x):
        total = 0
        for i in range(0, x):
            total = total + h[i]*(t**(x-(i+1)))
        return total

    def returnValue(self, z):
        AH = np.squeeze(np.asarray(np.linalg.solve(self.generatePolynomial(self.Xint), self.Yint)))
        return self.f(z, np.squeeze(np.asarray(AH)), AH.shape[0])

    def show(self, xlimit, ylimit, dots):
        self.graph(self.returnValue, "Polinomio", xlimit, ylimit, dots)

class Lagrange(InterpolationMethod):
    def __init__(self, x, y):
        InterpolationMethod.__init__(self, x,y)

    def returnValue(self, z):
        n = self.Xint.shape[0] #size
        L = 0 # value
        for i in range(0, n):
            num = 1
            den = 1
            for j in range(0,n):
                if(j != i):
                    num = num * (z-self.Xint[j])
                    den = den * (self.Xint[i] - self.Xint[j])
            L = L + self.Yint[i] * (num/den)
        return L

    def show(self, xlimt, ylimit, dots):
        self.graph(self.returnValue, "Lagrange", xlimt, ylimit, dots)

class newtonPolinomial(InterpolationMethod):
    def __init__(self, x,y):
        InterpolationMethod.__init__(self, x,y)

    def delta(self, x, y):
        i = x.shape[0]
        m = np.zeros(shape=(i,i))

        m[0:i, 0] = y[:]
        for a in range(1, i):
            for b in range(0, i-a):
                m[b,a] = round(m[b+1,a-1]-m[b, a-1], 8)/round(x[b+a]-x[b], 8)
        return m

    def returnValue(self, x):
        m = self.delta(self.Xint, self.Yint)
        result = self.Yint[0]
        for i in range(1, self.Yint.shape[0]):
            internal = m[0, i]
            mult:float = 1
            for j in range(0, i):
                mult = mult * (x-self.Xint[j])
            internal = internal *mult
            result = result + internal
        return result

    def show(self, xlimit, ylimit, dots):
        self.graph(self.returnValue, "Newton", xlimit, ylimit, dots)

#######GregNewton############
class GregNewtonPolynomial(InterpolationMethod):
    def __init__(self, x, y):
        InterpolationMethod.__init__(self, x, y)

    def gregNewtonPoli(self, y):
        i = y.shape[0]
        m = np.zeros(shape=(i, i))

        m[:, 0] = y[:]

        for a in range(1, i):
            for b in range(0, i-a):
                m[b,a] = round(m[b+1,a-1]-m[b, a-1], 8)
        return m

    def fac(self, i:int):
        result = 1;
        while i > 1:
            result = result * i
            i = i-1;
        return result

    def returnValue(self, x):
        m = self.gregNewtonPoli(self.Yint)
        result = self.Yint[0]
        for i in range(1, self.Yint.shape[0]):
            internal = round(m[0, i], 7)/self.fac(i)*1.0
            mult:float = 1
            u:float = round(x - self.Xint[0], 7)/round(self.Xint[1]-self.Xint[0], 7)
            for j in range(0, i):
                mult = mult * round((u - j), 7)
            internal = internal *mult
            result = result + internal
        return result

    def show(self, xlimit, ylimit, dots):
        self.graph(self.returnValue, "Newton_Gregory", xlimit, ylimit, dots)

####GregNewton#################
polinear = Polinear(xnp,ynp)
lagrange = Lagrange(xnp, ynp)
newton = newtonPolinomial(xnp, ynp)
GregNewton = GregNewtonPolynomial(xnp, ynp)

polinear.show(0, 0.8, 100)
lagrange.show(0, 0.8, 100)
newton.show(0, 0.8, 100)
GregNewton.show(0, 0.8, 100)

polinear.setInput(xpq, ypq)
lagrange.setInput(xpq, ypq)
newton.setInput(xpq, ypq)
GregNewton.setInput(xpq, ypq)

start_time = time.time()
polinear.show(1950, 1980, 100)
elapsed_time = time.time() - start_time

print("time Interpol: ")
print(elapsed_time)

start_time = time.time()
lagrange.show(1950, 1980, 100)
elapsed_time = time.time() - start_time

print("time: ")
print(elapsed_time)
newton.show(1950, 1980, 100)
GregNewton.show(1950, 1980, 100)

polinear.setInput(xp, yp)
lagrange.setInput(xp, yp)
newton.setInput(xp, yp)
GregNewton.setInput(xp, yp)

print(polinear.returnValue(1958))
print(lagrange.returnValue(1958))
print(newton.returnValue(1958))
print(GregNewton.returnValue(1958))

polinear.setInput(xpo, ypo)
lagrange.setInput(xpo, ypo)
newton.setInput(xpo, ypo)
GregNewton.setInput(xpo, ypo)

print(polinear.returnValue(1988))
print(lagrange.returnValue(1988))
print(newton.returnValue(1988))
print(GregNewton.returnValue(1988))

print("")

polinear.setInput(xpq, ypq)
lagrange.setInput(xpq, ypq)
newton.setInput(xpq, ypq)
GregNewton.setInput(xpq, ypq)

print(polinear.returnValue(1958))
print(lagrange.returnValue(1958))
print(newton.returnValue(1958))
print(GregNewton.returnValue(1958))

polinear.setInput(xpoq, ypoq)
lagrange.setInput(xpoq, ypoq)
newton.setInput(xpoq, ypoq)
GregNewton.setInput(xpoq, ypoq)

print(polinear.returnValue(1988))
print(lagrange.returnValue(1988))
print(newton.returnValue(1988))
print(GregNewton.returnValue(1988))

print()
