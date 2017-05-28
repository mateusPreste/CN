import matplotlib.pyplot as plt
import numpy as np

yp=np.array([139620, 173706, 311622])
xp=np.array([1950, 1960, 1970])

ypo=np.array([311622, 633392, 1011500])
xpo=np.array([1970, 1980, 1991])

ypq=np.array([139620, 173706, 311622, 633392])
xpq=np.array([1950, 1960, 1970, 1980])

ypoq=np.array([311622, 633392, 1011500, 1403793])
xpoq=np.array([1970, 1980, 1991, 2000])

def polinear(M):
    matx = np.zeros((M.shape[0], M.shape[0]))
    for i in range (0, M.shape[0]):
        for k in range(0, M.shape[0]):
            matx[i][k] = pow(M[i],M.shape[0]-(k+1))
    return matx

def f(t:int, h, x):
    total = 0

    for i in range(0, x):
        asf = x-(i+1)
        total = total + h[i]*(t**(x-(i+1)))
    return total

def interpol(z, x, y):
    AH = np.squeeze(np.asarray(np.linalg.solve(polinear(x), y)))
    return f(z, np.squeeze(np.asarray(AH)), AH.shape[0])

def lagrange(z, x, y):
    n = x.shape[0] #size
    L = 0 # value
    for i in range(0, n):
        num = 1
        den = 1
        for j in range(0,n):
            if(j != i):
                num = num * (z-x[j])
                den = den * (x[i] - x[j])
        L = L + y[i] * (num/den)
    return L

xnp = np.array([ 0.0, 0.2, 0.3, 0.5,  0.6])
ynp = np.array([1.008, 1.064, 1.125, 1.343, 1.512])

def newtonPoli(x, y):
    i = x.shape[0]
    m = np.zeros(shape=(i,i))

    m[0:i, 0] = y[:]

    for a in range(1, i):
        for b in range(0, i-a):
            m[b,a] = round(m[b+1,a-1]-m[b, a-1], 8)/round(x[b+a]-x[b], 8)
    return m

def returnPol(x, xp, yp):
    m = newtonPoli(xp, yp)
    result = yp[0]
    for i in range(1, yp.shape[0]):
        internal = m[0, i]
        mult:float = 1
        for j in range(0, i):
            mult = mult * (x-xp[j])
        internal = internal *mult
        result = result + internal
    return result

def PolinomioNewton(size, xnp, ynp):
    t = np.linspace(0, size, size)
    y = np.zeros(len(t))

    for i in range(len(t)):
        #y[i] = interpol(t[i], xp, yp)
        y[i] = returnPol(t[i], xnp, ynp)
    plt.plot(t, y)
    plt.show()

def showInterpol(xp, yp):
    t = np.linspace(0, 200, 200)
    y = np.zeros(len(t))

    for i in range(len(t)):
        y[i] = interpol(t[i], xp, yp)
    plt.plot(t, y)
    plt.show()

PolinomioNewton(200, xnp, ynp)
showInterpol(xnp, ynp)

print(interpol(1958, xp, yp))
print(lagrange(1958, xp, yp))
print(returnPol(1958, xp, yp))
print(interpol(1988, xpo, ypo))
print(lagrange(1988, xpo, ypo))
print(returnPol(1988, xpo, ypo))

print("")

print(interpol(1958, xpq, ypq))
print(lagrange(1958, xpq, ypq))
print(returnPol(1958, xpq, ypq))
print(interpol(1988, xpoq, ypoq))
print(lagrange(1988, xpoq, ypoq))
print(returnPol(1988, xpoq, ypoq))

print()
