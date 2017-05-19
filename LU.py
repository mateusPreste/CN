import numpy as np

def LU(U):
    L = np.zeros(U.shape)
    
    for i in range (0, U.shape[0]): # iterate vector
        pivot = U[i][i] # pivot in diagonal
        L[i][i] = 1 # fill main diagonal equal to 1
        
        for j in range (i+1, U.shape[0]):
            m = U[j][i]/pivot #find multi
            L[j][i] = m
            U[j][:] = U[j][:] - (m*U[i][:])
            
    y = U[:,U.shape[1]-1] 
    
    return supFind( np.delete(U, U.shape[1]-1, 1), y)

def infFind(A, result): #A -> matrix and result -> result vector
    coef = np.zeros(result.shape, dtype=np.float64) #coeficients vector
    for i in range (0, A.shape[0]):
        sum = np.float64(0.0)
        for k in range (0, i):
            sum = sum + (A[i][k]*coef[k])
        
        coef[i] = (result[i] - sum)/(A[i][i])
    return coef

def supFind(A, result): #A -> matrix and result -> result vector
    coef = np.zeros(result.shape, dtype=np.float64) #coeficients vector
    for i in range (A.shape[0]-1, -1, -1):
        sum = np.float64(0.0)
        for k in range (A.shape[1]-1, i, -1):
            sum = sum + (A[i][k]*coef[k])
        
        coef[i] = (result[i] - sum)/(A[i][i])
    return coef

uu = np.array( [ [1.,6,2,4,8],
                [3,19,4,15,25],
                [1,4,8,-12,18],
                [5,33,9,3,72] ])

print(LU(uu))
