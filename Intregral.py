import numpy as np
import matplotlib.pyplot as plt

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
        plt.ylim([0, self.Yint(np.argmax(self.Yint))*1.4])
        plt.show()

class Polinear(InterpolationMethod):
    def __init__(self, x, y):
        InterpolationMethod.__init__(self,x,y)

    def pow(self, a:int, b):
        result = 1
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

#Parent Method
class IntegrationMethod():
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

def f1sobreX(x):
    return 1/x;

class Trapezio(IntegrationMethod):
    def __init__(self, x, y):
        IntegrationMethod.__init__(self,x,y)

    def f(self, fun, a, b, n): #r é o numero de repetições
        h = (b-a)
        soma = fun(a)
        for i in range(1, n-1):
            soma += 2*(fun(a+(h*i)))
        soma += fun(b)
        return (h/2)*(soma)

    def returnValue(self, z):
        return self.f()

    def show(self, xlimit, ylimit, dots):
        self.graph(self.returnValue, "Polinomio", xlimit, ylimit, dots)
        
global_Xint = None
global_Yint = None
    
def trapezio(fun, a, b, n): #r é o numero de repetições
    global global_Xint
    global global_Yint
    global_Xint = np.zeros(n+1)
    global_Yint = np.zeros(n+1)
    global_Xint[0] = a
    global_Xint[n] = b
    global_Yint[0] = fun(a)
    global_Yint[n] = fun(b)
    h = (b-a)/n
    soma = (fun(a)+fun(b))
    for i in range(1, n):
        soma += 2*(fun(a+ (i*h)))
        global_Xint[i] = a+(i*h)
        global_Yint[i] = fun(global_Xint[i])
    return (h/2)*(soma)
    
def umTerco(fun, a, b, n): #r é o numero de repetições
    global global_Xint
    global global_Yint
    global_Xint = np.zeros(n+1)
    global_Yint = np.zeros(n+1)
    global_Xint[0] = a
    global_Xint[n] = b
    global_Yint[0] = fun(a)
    global_Yint[n] = fun(b)
    h = (b-a)/n
    soma = fun(a)
    for i in range(1, n):
        if(i%2 == 0):
            soma += 2*(fun(a+(h*i)))
        else:
            soma += 4*(fun(a+(h*i)))
        global_Xint[i] = a+(i*h)
        global_Yint[i] = fun(global_Xint[i])
    soma += fun(b)
    return (h/3)*(soma)
    
def tresOitavos(fun, a, b, n): #r é o numero de repetições
    h = (b-a)/n
    soma = fun(a)
    for i in range(1, n-1):
        if(i%3 == 0):
            soma += 2*(fun(a+(h*i)))
        else:
            soma += 3*(fun(a+(h*i)))
    soma += fun(b)
    return ((3*h)/8)*(soma)

def graph(function, legend, xlimit, ylimit, dots):
        t = np.linspace(xlimit, ylimit, dots)
        y = np.zeros(len(t))

        for i in range(len(t)):
            y[i] = function(t[i])
        
        plt.plot(t, y, label=legend)
        plt.plot(global_Xint, global_Yint, 'bo')
        plt.plot(global_Xint, global_Yint, 'ro', label = 'Pontos inseridos')
        plt.plot(global_Xint, global_Yint, '-o')
        plt.title(legend)
        plt.legend();
        #plt.ylim([0, global_Yint[Yint.shape[0]-1]*1.4])
        plt.show()
        
for i in range(1, 31):
    print("Trapezio: ")
    print(trapezio(f1sobreX, 1, 7, i))
    graph(f1sobreX, "Trapezio", 1, 7, 50)

print()
print(trapezio(f1sobreX, 1, 7, 2))
graph(f1sobreX, "Trapezio", 1, 7, 50)
print()
print(trapezio(f1sobreX, 1, 7, 30))
graph(f1sobreX, "Trapezio", 1, 7, 50)
print()
#print(trapezio(f1sobreX, 1, 7, 30000))
#graph(f1sobreX, "Trapezio", 1, 7, 50)
print()




print("1/3: ")
print(umTerco(f1sobreX, 1, 7, 1))
print(global_Xint)
print(global_Yint)


print()
print(umTerco(f1sobreX, 1, 7, 2))
polinom = Polinear(global_Xint, global_Yint)
polinom.show(1, 7, 200)

print()
print(umTerco(f1sobreX, 1, 7, 30))

print()
#print(umTerco(f1sobreX, 1, 7, 30000))

print("3/8: ")
print(tresOitavos(f1sobreX, 1, 7, 1))
print()
print(tresOitavos(f1sobreX, 1, 7, 2))
print()
print(tresOitavos(f1sobreX, 1, 7, 30))

print()
#print(tresOitavos(f1sobreX, 1, 7, 30000))
