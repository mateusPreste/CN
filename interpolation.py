from matplotlib.pyplot import *
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

def f(t, h, x):
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
    
print(interpol(1958, xp, yp))
print(lagrange(1958, xp, yp))
print(interpol(1988, xpo, ypo))
print(lagrange(1988, xpo, ypo))

print("")

print(interpol(1958, xpq, ypq))
print(lagrange(1958, xpq, ypq))
print(interpol(1988, xpoq, ypoq))
print(lagrange(1988, xpoq, ypoq))

    
t = np.linspace(1950, 1970, 2000)
y = np.zeros(len(t))

for i in range(len(t)):
    y[i] = interpol(t[i], xp, yp)
plot(t, y)


show()
