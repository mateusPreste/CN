import matplotlib.pyplot as plt
import numpy as np
import time
import math as mt

y = np.array([139620, 173706, 311622, 633392, 1011500, 1403796])
x = np.array([1950, 1960, 1970, 1980, 1991, 2000])

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

    def pow(self, a:int, b):
        result = 1
        result = np.int64(result)
        for i in range(0, b):
            result = result * a
        return result

    def Matrix(self, M):
        matx = np.zeros((M.shape[0], M.shape[0]))
        for i in range (0, M.shape[0]):
            for k in range(0, M.shape[0]):
                matx[i][k] = self.pow(M[i],M.shape[0]-(k+1))
        return matx

    def f(self, t, h, x):
        total = 0
        for i in range(0, x):
            total = total + h[i]* ( t** (x-(i+1)) )
        return total

    def returnValue(self, z):
        AH = np.squeeze(np.asarray(np.linalg.solve(self.Matrix(self.Xint), self.Yint)))

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
questao1X = np.array([1250, 1000, 750])
questao1Y = np.array([25, 15, 10])



polinear = Polinear(questao1X,questao1Y)
lagrange = Lagrange(questao1X,questao1Y)
newton = newtonPolinomial(questao1X,questao1Y)
GregNewton = GregNewtonPolynomial(questao1X,questao1Y)


#Questão 1

def questao(value):
    polinear.returnValue(value)
    start_time = time.time()
    resultP = polinear.returnValue(value)
    elapsed_time = time.time() - start_time

    lagrange.returnValue(value)
    start_timeL = time.time()
    resultL = lagrange.returnValue(value)
    elapsed_timeL = time.time() - start_timeL

    newton.returnValue(value)
    start_timeN = time.time()
    resultN = newton.returnValue(value)
    elapsed_timeN = time.time() - start_timeN

    GregNewton.returnValue(value)
    start_timeG = time.time()
    resultG = GregNewton.returnValue(value)
    elapsed_timeG = time.time() - start_timeG
    print("polinear   "+ str(resultP) + " tempo: "+ str(elapsed_time))
    print("lagrange   "+ str(resultL) + " tempo: "+ str(elapsed_timeL))
    print("newton     "+ str(resultN) + " tempo: "+ str(elapsed_timeN))
    print("GregNewton "+ str(resultG) + " tempo: "+ str(elapsed_timeG))

def setRange(polinear, lagrange, newton, GregNewton, x, y, a,b):
    polinear.setInput(x[a:b],y[a:b])
    lagrange.setInput(x[a:b],y[a:b])
    newton.setInput(x[a:b],y[a:b])
    GregNewton.setInput(x[a:b],y[a:b])

y = np.array([139620, 173706, 311622, 633392, 1011500, 1403796])
x = np.array([1950, 1960, 1970, 1980, 1991, 2000])

print("Questao 1")
print()
print("850 metros")
questao(850)
print()
print("1050 metros")
questao(1050)
print()
print("1200 metros")
questao(1200)

print()
print("Questao 2")
print()
setRange(polinear, lagrange, newton, GregNewton, x, y, 0, 3)
print("1958 3 pontos")
questao(1958)

setRange(polinear, lagrange, newton, GregNewton, x, y, 2, 5)

print()
print("1988 3 pontos")
questao(1988)
print()

setRange(polinear, lagrange, newton, GregNewton, x, y, 0, 4)
print("1958 4 pontos")
questao(1958)

setRange(polinear, lagrange, newton, GregNewton, x, y, 2, 6)

print()
print("1988 4 pontos")
questao(1988)
print()

def function3(x):
    return 10*x**4+2*x+1

def function4(x):
    return mt.sin(x)+2*x+1

x = np.array([0.1, 0.2, 0.3])
y = np.array([function3(0.1), function3(0.2), function3(0.3)])
x4 = np.array([0.1, 0.2, 0.3, 0.4])
y4 = np.array([function4(0.1), function4(0.2), function4(0.3), function4(0.4)])

setRange(polinear, lagrange, newton, GregNewton, x, y, 1, 3)

print("Questao 3")
print()
print("P2(0.15)")
questao(0.15)
print("erro absoluto")
print(abs(polinear.returnValue(0.15)-function3(0.15)))
print("erro relativo")
print(str(abs(polinear.returnValue(0.15)-function3(0.15))/function3(0.15)*100)+" %")

print()
print("Questao 4")
print()
setRange(polinear, lagrange, newton, GregNewton, x4, y4, 1, 3)
print("L2(0.15)")
questao(0.15)
print("erro absoluto")
print(abs(lagrange.returnValue(0.15)-function4(0.15)))
print("erro relativo")
print(str(abs(lagrange.returnValue(0.15)-function4(0.15))/function4(0.15)*100)+" %")

print()

setRange(polinear, lagrange, newton, GregNewton, x4, y4, 1, 4)

print("L3(0.15)")
questao(0.15)
print(function4(0.15))
print("erro absoluto")
print(abs(lagrange.returnValue(0.15)-function4(0.15)))
print("erro relativo")
print(str(abs(lagrange.returnValue(0.15)-function4(0.15))/function4(0.15)*100)+" %")

'''
polinear.show(0, 0.8, 100)
lagrange.show(0, 0.8, 100)
newton.show(0, 0.8, 100)
GregNewton.show(0, 0.8, 100)

polinear.setInput(x[0:4], y[0:4])
lagrange.setInput(x[0:4], y[0:4])
newton.setInput(x[0:4], y[0:4])
GregNewton.setInput(x[0:4], y[0:4])

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












#declarando x e y
y = np.array([139620, 173706, 311622])
x = np.array([1950, 1960, 1970])

#declarando
polinear = Polinear(xnp,ynp)
lagrange = Lagrange(xnp, ynp)
newton = newtonPolinomial(xnp, ynp)
GregNewton = GregNewtonPolynomial(xnp, ynp)

polinear.setInput(x, y) #usando x e y em
lagrange.setInput(x, y) #usando x e y em
newton.setInput(x, y) #usando x e y em
GregNewton.setInput(x, y) #usando x e y em

print(polinear.returnValue(1988)) #interpolação de 1988 usando metodo polinear
print(lagrange.returnValue(1988)) #interpolação de 1988 usando metodo lagrange
print(newton.returnValue(1988)) #interpolação de 1988 usando metodo newton
print(GregNewton.returnValue(1988)) #interpolação de 1988 usando metodo greg e newton



'''





print()
