#Random Number Generator
class pRNG():
    '''My class to store seed value and different types of pRNG'''
    seed = 101
    def LCG(low=0,high=1):
        '''Gives a random number between low and high'''
        a = 1103515245
        c = 12345
        m = 32768
        x1 = ((a*(pRNG.seed)+c)%m)
        pRNG.seed = x1 #changing the seed value every time the rng is called
        val = low + abs(high-low)*x1/m
        return val

##Linear Algebra

from math import *

# function to print a matrix aesthetically
def Print_Matrix(a):[print((' \t').join(str(k) for k in a[i]),'\n') for i in range(len(a))]
    
#Function to read the arrays from an external .txt file
def read_arr(path):
    contents = {}
    with open(path,'r') as f:
        data = f.readlines()
    #print(data)
    def convert_frac(s):
        if '/' in s:
            return float(s.split('/')[0])/float(s.split('/')[1])
        else:
            return float(s)
    for line in data:
        lin_lst = line.strip().split()
        #print(lin_lst)
        if lin_lst[0] == '#':
            name = lin_lst[1]
            contents[name] = []
        else:
            if 'matrix' in name:
                contents[name].append([convert_frac(i) for i in lin_lst])
            elif 'vector' in name:
                contents[name] = [convert_frac(i) for i in lin_lst]
    return contents

#Function to read x and y values from a file
def read_file(filename: str,delimiter: str = '\t'):

    matrices = []
    current_matrix = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()  
            if not line or line.startswith("#"):
                if current_matrix: 
                    matrices.append(current_matrix)
                    current_matrix = []  
                continue
            
            try:
                row = [float(num) for num in line.split(delimiter)]
                current_matrix.append(row)
            except ValueError:
                # print("Skipping non-numeric line:", line)
                pass
        if current_matrix:
            matrices.append(current_matrix)
    return matrices

#let's say I have a matrix A and a matrix B:
# then their augmented matrix code is:
def aug_mtrx(a,b):
    n = len(a) #number of rows of matrix A
    if type(b[0])!=list:
        aug = [a[i]+[b[i]] for i in range(n)]
    else:
        aug = [a[i]+b[i] for i in range(n)]
    return aug

#defining swapping function in order to create reduce row echelon form
def swap(a,i):
    swap = 0
    if abs(a[i][i])<1e-3:
        for j in range(i+1,len(a)-1):
            if not (abs(a[i][j])<1e-3 and abs(a[j][j])<1e-3):
                swap += 1
                a[i],a[j]=a[j],a[i]
        if swap == 0:
            return None
    return swap

#Gauss-Jordan method

#Now the code for Gauss-Jordan method:
def Gauss_Jordan_Solve(a,b,r=4):
    n = len(a)
    A = aug_mtrx(a,b)
    for i in range(n):
        if swap(A,i)==None:
            print("solution is not possible")
            return None
        A[i] = [x/A[i][i] for x in A[i]]
        for j in range(n):
            if j!=i:
                A[j] = [A[j][k]-A[j][i]*A[i][k] for k in range(n+1)]
    A = [round(ele[-1],r) for ele in A]
    return A

#Code to form an Identity matrix of order n:
def Identity(n):
    I_n = [n*[0] for i in range(n)]
    for i in range(n):
        I_n[i][i] = 1
    return I_n

##Code to find the inverse of a matrix using Gauss-Jordan method
def LU_Inverse(a,r_factor = 5):
    n = len(a)
    b = Identity(n)
    A = aug_mtrx(a,b)
    #Print_Matrix(A)
    for i in range(n):
        if swap(A,i)==None:
            print("solution is not possible")
            return None
        A[i] = [x/A[i][i] for x in A[i]]
        for j in range(n):
            if j!=i:
                A[j] = [A[j][k]-A[j][i]*A[i][k] for k in range(2*n)]
    #Print_Matrix(A)
    inverse = [[round(x,r_factor) for x in ele] for ele in (A[i][-n:] for i in range(n))]
    #Print_Matrix(inverse)
    return inverse

#Code to find the determinant of matrix a:
def Determinant(A,r=4):
    n = len(A)
    a = [A[i][:] for i in range(n)]
    count = 0
    for i in range(n):
        swap_i = swap(a,i)
        if swap_i==None:
            print("Solution is not possible!")
            return None
        count += swap_i
        for j in range(i+1,n):
            a[j] = [a[j][k]-a[j][i]*(a[i][k]/a[i][i] )for k in range(n)]
    #Print_Matrix(a)
    det = (-1)**count
    for l in range(n):
        det *= a[l][l]
    return round(det,r)

#L-U Decomposition

#Function to decompose a square matrix into lower and upper matrix
def LU_decompose(A):
    n = len(A)
    U = [A[i][:] for i in range(n)]
    L = Identity(n)
    for i in range(n):
        for j in range(i+1,n):
            L[j][i] = U[j][i]/U[i][i]
            U[j] = [U[j][k]-L[j][i]*U[i][k] for k in range(n)]
    return L,U

#Function to perform forward substitution
def forward_sub(L,B):
    n = len(L)
    b = [B[i][0] if type(B[i])==list else B[i] for i in range(n)]
    #print(b)
    y = n*[0]
    y[0] = b[0]
    for i in range(1,n):
        y[i] = b[i] - sum(L[i][j]*y[j] for j in range(i))
    return y

#Function to perform backward substitution
def backward_sub(U,y):
    n = len(U)
    x = n*[0]
    x[-1] = y[-1]/U[-1][-1]
    for i in range(n-2,-1,-1):
        x[i] = (y[i] - sum(U[i][j]*x[j] for j in range(i+1,n)))/U[i][i]
    return x

#Function to solve a system of linear equation using LU-Decomposition
def LU_solve(a,b,r=6):
    L,U = LU_decompose(a)
    y = forward_sub(L,b)
    x = backward_sub(U,y)
    return [round(i,r) for i in x]

#Cholesky Decomposition

#Function to check whether the given matrix is symmetric or not
def isSymmetric(a):
    n = len(a)
    for i in range(n):
        for j in range(n):
            if a[i][j] != a[j][i]:
                return False
    return True
    
#Function to take a transpose of given square matrix
def transpose(a):
    n = len(a)
    tp = [n*[0] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            tp[i][j] = a[j][i]
    return tp

#Cholesky decomposition function returning an Lower and an Upper matrix
def Cholesky_decompose(A):
    if isSymmetric(A) == False:
        print('Matrix is not symmetric')
        return None
    n = len(A)
    L = [n*[0] for _ in range(n)]
    U = [n*[0] for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            cal = A[i][j] - sum(L[i][k]*L[j][k] for k in range(j))
            if j == i:  
                L[i][j] = sqrt(cal)
            else:
                L[i][j] = cal/L[j][j]
            U[j][i] = L[i][j]
    return L,U

#Function to perform forward substitution for cholesky decomposition
def forward_sub1(L,B):
    n = len(L)
    b = [B[i][0] if type(B[i])==list else B[i] for i in range(n)]
    #print(b)
    y = n*[0]
    for i in range(n):
        y[i] = (b[i] - sum(L[i][j]*y[j] for j in range(i)))/L[i][i]
    return y

#Function to solve system f equations using cholesky decomposition
def Cholesky_solve(A,B,r=6):
    L,U = Cholesky_decompose(A)
    Y = forward_sub1(L,B)
    X = backward_sub(U,Y)
    return [round(ele,r) for ele in X]

#Jacobi Method

#Function to check whether the given matrix is diagonally dominant
def isDiaDominant(A):
    for i in range(len(A)):
        if not abs(A[i][i]) >= sum(abs(A[i][j]) for j in range(len(A)) if j!=i):
            return False
    return True

# make a matrix diagonally dominant 
def diagDM (A,B):
    n = len(A)
    if isDiaDominant(A) != True:
        for i in range(n):
            add = sum(abs(A[i][j]) for j in range(n) if j!=i)
            if abs(A[i][i]) < add:
                for k in range(i,n):
                    A[i],A[k] = A[k],A[i]
                    B[i],B[k] = B[k],B[i]
                    s = sum(abs(A[i][m]) for m in range(n) if m!=i)
                    if A[i][i] > s: break
    return A , B

#Function to find solution of system of linear equations using Jacobi Method
def Jacobi_method(A,B,xj,e=6):
    n = len(A)
    xi = n*[0]
    condition = True
    count = 0
    while condition:
        for i in range(n):
            xi[i] = round((1/A[i][i])*(B[i]-sum(A[i][j]*xj[j] for j in range(n) if j!=i)),6)
        condition = False not in [abs(xi[k]-xj[k])>10**(-1*e) for k in range(n)]
        xj = xi[:]
        count += 1
        #print(f"Iteration number: {count}, {[f'{i:{1}.{e}}' for i in xi]}\n")
    xi = [round(i,e) for i in xi]
    print(f"No of Iterations: {count}\n")
    return xi

#Gauss-Seidel Method
#Function to find solution of system of linear equations using Gauss-Seidel Method
def Gauss_Seidel(A,B,xi,e=6):
    #Finding length of a    
    n = len(A)
    condition = True
    count = 0
    while condition:
        xj = xi[:]
        for i in range(n):
            xi[i] = (B[i]-sum(A[i][j]*xi[j] for j in range(i))-sum(A[i][j]*xj[j] for j in range(i+1,n)))/A[i][i]
        count += 1
        condition = sum([abs(xi[k]-xj[k])>10**(-1*e) for k in range(n)])==n
        #print(f"Iteration number: {count}, {[f'{i:{1}.{e}}' for i in xi]}\n")
    xi = [round(i,e) for i in xi]
    print(f"No of Iterations: {count}\n")
    return xi

def cgrad(A, x0, b, max_iterations, tolerance):
    x = x0
    r = b - A.dot(x)
    p = r
    rsold = r.dot(r)
    from numpy import dot
    for i in range(max_iterations):
        Ap = A.dot(p)
        alpha = rsold / p.dot(Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = r.dot(r)
        
        if np.sqrt(rsnew) < tolerance:
            break
        
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    
    return x

def cg_fly( x0, b, max_iterations, tolerance):
    def delta(x, y):
        if x == y:
            return 1
        else:
            return 0
    def mat(x, y):
        return (delta(x+1, y) + delta(x-1, y) - 2*delta(x, y))*0.5 + 0.04* delta(x, y)
        # return (delta(x+1, y) + delta(x-1, y) + 2*delta(x, y))*0.5
        # return delta(x,y)
    def dot2(x):
        n = len(x)
        r = np.zeros(n)
        for row in range(n):
            for i in range(n):
                    r[row] += mat(row,i)*x[i]
        return np.array(r)
    
    x = x0
    r = b - dot2(x)
    p = r
    rsold = r.dot(r)
    res = []

    for i in range(max_iterations):
        Ap = dot2(p)
        alpha = rsold / p.dot(Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = r.dot(r)
        if np.sqrt(rsnew) < tolerance:
            break
        res.append(rsold)
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    return x, np.array(res)


##Functions for root finding

#Function to find root using Bisection method
def Bisection_root(fn,interval,e=6):
    a,b = interval
    if (fn(a) * fn(b) >= 0):
        print("You have not assumed right interval as the sign of value of function at these points is not opposite\n")
        return None
    c = a
    n = 0
    while (b-a) >= 10**-e and abs(fn(c)) > 10**-e:
        n += 1
        c = (a+b)/2
        #if abs(fn(c)) <= 10**-e:
            #break
        if (fn(c)*fn(a)<0):
            b = c
        else:
            a = c
        print(f"Iteration no: {n} \troot -> {c:{1}.{e}f}")
    print(f"\nThe root in the given interval converges to {c:{1}.{e}f} and the value of function is {fn(c):{1}.{e}f}")
    print("No. of iterations = ",n)
    #return c

#Function to show bracketing
def Bracketing(fn,interval,k=1,max_iter = 12):
    a,b = interval
    n = 1
    while not fn(a)*fn(b) < 0 and n <= max_iter:
        if abs(fn(a)) < abs(fn(b)):
            a -= k*(b-a)
        elif abs(fn(a)) > abs(fn(b)):
            b += k*(b-a)
        print(f"N: {n} \t\tBracketing: \t({a}, {b})")
        n += 1
    print("\nTotal no of iterations = ",n-1)
    return [a,b]

#Function to find root using Regula-Falsi method
def Regula_falsi(fn,interval,e=6):
    a,b = interval
    if (fn(a) * fn(b) >= 0):
        print("You have not assumed right interval\n")
        return None
    n = 0
    c=a
    c1 = c-1
    while abs(c-c1) > 10**-e or abs(fn(c)) > 10**-e:
        n += 1
        c = (a * fn(b) - b * fn(a))/ (fn(b) - fn(a))
        #if abs(fn(c)) <= 10**-e:
            #break
        if (fn(c)*fn(a)<0):
            b = c
        elif (fn(c)*fn(b)<0):
            a = c
        c1 = c
        print(f"Iteration no: {n} \troot -> {c:{1}.{e}f}")
    print(f"\nThe root in the given interval converges to {c:{1}.{e}f} and the value of function is {fn(c):{1}.{e}f}")
    print("Total no of iterations = ",n)
    #return c

#Function to find root using Newton_raphson method
def Newton_Raphson(fn,d_fn,x = 0.5,e=6):
    '''fn: the function of which we want to find the root,
       d_fn: the derivative of the function
       x: initial guess for the root'''
    h = fn(x)/d_fn(x)
    n = 0
    while abs(h)>10**-e:
        x -= h
        h = fn(x)/d_fn(x)
        n += 1
        #print(f"Iteration no: {n} \troot -> {x:{1}.{e}f}")
    #print(f"The root converges to {x:{1}.{e}f} and the value of function is {fn(x):{1}.{e}f}")
    print("\nTotal no of iterations = ",n)
    return round(x,e)

#Laguerre's Method

# Function to find the value of a polynomial at a given 'x'
# using its list of coefficiewnts
def Pn(x,c_lst):
    n = len(c_lst)-1
    add = 0
    for c in c_lst:
        add += c*(x**n)
        n -= 1
    return add

# Function to get a list of coeff of derivative of a polinomial
def derivative(c_lst):
    n = len(c_lst)-1
    return [c_lst[i]*(n-i) for i in range(n)]

# to get a list of coeff after deflating a polinomial
def deflation(root,c_lst):
    coeff = c_lst[:-1]
    k = 0
    for i in range(len(coeff)):
        coeff[i] += k
        k = root*coeff[i]
    return coeff

# Laguerre method to find all the roots of a polinomial
def Laguerre(coeff,r_lst=[],guess=1,e=4):
    c_lst = coeff[:]
    b0 = guess
    n = len(c_lst)-1
    for r in r_lst:
        c_lst = deflation(r,c_lst)
    while len(r_lst) != n:
        if abs(Pn(b0,c_lst)) < 10**-e:
            r_lst.append(round(b0,e))
            c_lst = deflation(b0,c_lst)
            if len(c_lst) == 2:
                r_lst.append(round(-c_lst[1],e))
                break
            continue
        while not abs(Pn(b0,c_lst)) < 10**-e:
            c1 = derivative(c_lst)
            c2 = derivative(c1)
            G = Pn(b0,c1)/Pn(b0,c_lst)
            H = G*G - Pn(b0,c2)/Pn(b0,c_lst)
            a = n/(G+abs(((n-1)*(n*H-G*G))**0.5))
            b1 = b0 - a
            if abs(b1-b0) < 10**-e:
                b0 = b1
                break
            b0 = b1
    return r_lst

##Data Interpolation and Fitting

#Lagrange Interpolation
def Lagrange(x,X_lst,Y_lst,N):
    y = 0
    for i in range(N+1):
        prod = Y_lst[i]
        for k in range(N+1):
            if i!=k:
                prod *= (x-X_lst[k])/(X_lst[i]-X_lst[k])
        y += prod
    return y

#Least Square Fitting

# for linear relation
def Linear_regression(X_lst,Y_lst,sigma=[]):
    N = len(X_lst)
    if len(sigma) == 0: sigma = N*[1]
    Sy = sum(Y_lst[i]/sigma[i]**2 for i in range(N))
    S = sum(1/sigma[i]**2 for i in range(N))
    Sx = sum(X_lst[i]/sigma[i]**2 for i in range(N))
    Sxy = sum(X_lst[i]*Y_lst[i]/sigma[i]**2 for i in range(N))
    Sxx = sum((X_lst[i]/sigma[i])**2 for i in range(N))
    Syy = sum((Y_lst[i]/sigma[i])**2 for i in range(N))
    delta = S*Sxx - Sx*Sx
    c = (Sxx*Sy - Sx*Sxy)/delta
    m = (Sxy*S - Sx*Sy)/delta
    var_c = Sxx/delta
    var_m = S/delta
    pearson_r = (Sxy*Sxy)/(Sxx*Syy)
    return m,c,var_m,var_c,pearson_r

#for polynomial of degree k
def LeastSquare_polyfit(X,Y,k):
    n = len(X)
    A = [(k+1)*[0] for _ in range(k+1)]
    B = (k+1)*[0]
    for i in range(k+1):
        B[i] = sum((X[l]**i)*Y[l] for l in range(n))
        A[0][i] = sum(X[l]**i for l in range(n))
    for j in range(1,k+1):
        A[j] = A[j-1][1:] + [sum((X[l]**(k+j)) for l in range(n))]
    solution = Gauss_Jordan_Solve(A,B)
    return solution

##Numerucal Integration

#Function to find integration using Midpoint method
def midpoint_int(fn,a,b,n):
    h = abs(a-b)/n
    add = 0
    xi = a+h/2
    for i in range(n):
        add += h*fn(xi)
        xi += h
    return add

#Function to find integration using Trapezoidal method
def trap_int(fn,a,b,n):
    h = abs(a-b)/n
    add = 0
    xi = a
    for i in range(n+1):
        if i == 0 or i == n:
            add += (h/2)*fn(xi)
        else:
            add += h*fn(xi)
        xi += h
    return add

#Function to find integration using Simpson method
def Simpson_int(fn,x0,x2,n):
    h = abs(x2-x0)/n
    add = 0
    xi = x0
    for i in range(n+1):
        if i == 0 or i == n:
            add += (h/3)*fn(xi)
        elif i%2 == 0:
            add += (2*h/3)*fn(xi)
        else:
            add += (4*h/3)*fn(xi)
        xi += h
    return add

#Function to find integration using Monte-Carlo method
def Monte_Carlo(fn,a,b,N):
    from mylib import pRNG
    f_x = [fn(pRNG.LCG(a,b)) for _ in range(N)]
    avg = sum(f_x)/N
    s = (b-a)*avg
    var = (1/N)*sum(i*i for i in f_x) - (avg)**2
    return s, var

#RK4
import matplotlib.pyplot as plt

def RK4(xt,f_xy,x0,y0,h=1e-3):
    n = (xt-x0)//h
    X_arr = [x0]
    Y_arr = [y0]
    for i in range(n):
        k1 = h*f_xy(y0,x0)
        k2 = h*f_xy(y0 + k1/2,x0 + h/2)
        k3 = h*f_xy(y0 + k2/2,x0 + h/2)
        k4 = h*f_xy(y0 + k3,x0 + h)
        yt = y0 + (1/6)*(k1+2*k2+2*k3+k4)
        y0 = yt
        x0 += h
        X_arr.append(x0)
        Y_arr.append(y0)
    plt.plot(X_arr,Y_arr,'g-',label = 'final plot')
    plt.grid()
    plt.legend()
    plt.show()
    return X_arr, Y_arr

def RK4_Coupled_XY(d2ydx2, dydx, x0, y0, z0, xf, st):
    x = [x0]
    y = [y0]
    z = [z0]      # dy/dx
    n = int((xf-x0)/st)     # no. of steps
    for i in range(n-1):
        x.append(x[i] + st)
        k1 = st * dydx(x[i], y[i], z[i])
        l1 = st * d2ydx2(x[i], y[i], z[i])
        k2 = st * dydx(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        l2 = st * d2ydx2(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        k3 = st * dydx(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        l3 = st * d2ydx2(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        k4 = st * dydx(x[i] + st, y[i] + k3, z[i] + l3)
        l4 = st * d2ydx2(x[i] + st, y[i] + k3, z[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (l1 + 2*l2 + 2*l3 + l4)/6)
        
    plt.plot(x,y,'r-', label = 'Y vs T')
    plt.plot(x,z,'g-', label = 'V vs T')
    plt.grid()
    plt.legend()
    plt.show()
    #print('Line here is not a fitting of the polynomial. Has been added to aid the eye to track the points.')

    return x, y, z

def RK4_Coupled_XYZ(dxdt, dydt, dzdt, x0, y0, z0, t0, tf, st):
    x = [x0]
    y = [y0]
    z = [z0]
    t = [t0] 

    n = int((tf-t0)/st)     # no. of steps
    for i in range(n):
        t.append(t[i] + st)
        k1x = st * dxdt(x[i], y[i], z[i],t[i])
        k1y = st * dydt(x[i], y[i], z[i],t[i])
        k1z = st * dzdt(x[i], y[i], z[i],t[i])
        
        k2x = st * dxdt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        k2y = st * dydt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        k2z = st * dzdt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        
        k3x = st * dxdt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        k3y = st * dydt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        k3z = st * dzdt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        
        k4x = st * dxdt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        k4y = st * dydt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        k4z = st * dzdt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        

        x.append(x[i] + (k1x + 2*k2x + 2*k3x + k4x)/6)
        y.append(y[i] + (k1y + 2*k2y + 2*k3y + k4y)/6)
        z.append(z[i] + (k1z + 2*k2z + 2*k3z + k4z)/6)

    sphere = plt.axes(projection='3d')
    sphere.plot(x, y, z,'o-')
    
    return x, y, z, t

def RK4_boundary_value(y0, z0,x0, N, end, interpol):
    y_i = y0
    z_i = z0
    step = 0
    yl = [y0]
    zl = [z0]
    xl = [x0]
    h = (end-0)/N
    while step <= end:
        k1y = h*boundary_RK4_dyx(y_i, z_i, step)
        k1z = h*boundary_RK4_dzx(y_i, z_i, step)

        k2y = h*boundary_RK4_dyx(y_i + k1y/2, z_i+k1z/2, step+h/2)
        k2z = h*boundary_RK4_dzx(y_i + k1y/2, z_i+k1z/2, step+h/2)

        k3y = h*boundary_RK4_dyx(y_i+k2y/2, z_i+k2z/2, step+h/2)
        k3z = h*boundary_RK4_dzx(y_i+k2y/2, z_i+k2z/2, step+h/2)

        k4y = h*boundary_RK4_dyx(y_i+k3y, z_i+k3z, step+h)
        k4z = h*boundary_RK4_dzx(y_i+k3y, z_i+k3z, step+h)

        y_i += (k1y+2*k2y+2*k3y+k4y)/6
        z_i += (k1z+2*k2z+2*k3z+k4z)/6
        yl.append(y_i)
        zl.append(z_i)
        step += h
        xl.append(step)

    if interpol == 0:
        return (zl, yl, xl)
    else:
        return (z_i, y_i)

def Shooting_method(y_0, y_L, y_x, N):
    z1 = -0.5 #inital guess
    z1h, y1 = RK4_boundary_value(y_0, z1,0, y_x, 10, 1) #initialiing
    z2 = -2 #inital guess
    z2h, y2 = RK4_boundary_value(y_0, z2,0, y_x, 10, 1) #initialiing
    iter = 0
    while abs(y1-y_L) >= 0.001 and abs(y2-y_L)>=0.001 and iter <= 30:
        iter +=1
        znew = z2 + ((z1 - z2)*(y_L - y2))/(y1 - y2)
        znew2, ynew = RK4_boundary_value(y_0, znew,0, N, 10, 1)
        #print(ynew, znew)
        if abs(ynew - y_L) < 0.001:
            z, y, x = RK4_boundary_value(y_0, znew,0, N, 10, 0)
            break
        else:
            if ynew < y_L:
                z2 = znew
                y2 = ynew
            else:
                z1 = znew
                y1 = ynew
    plt.plot(x, y,'r',label="T vs position")
    for i in range(0, len(y)):
        if abs(y[i] - 100) < 0.1:#to get value of position at temperature 100
            out = x[i]
            break
    return out

def boundary_RK4_dyx(y,z,x):
    f = z
    return f

def boundary_RK4_dzx(y,z,x):
    f = 0.01*(y-20)
    return f

import numpy as np
import toolsar, math

def PDE_Solve_1d_heat(lx,Nx,lt,Nt,lower_x,tot_steps):
    step_arr = [10, 20, 50, 100, 200, 500, 1000, tot_steps]
    hx=(lx/Nx)
    ht=(lt/Nt)
    alpha=ht/(hx)**2
    V0=np.zeros(Nx+1)
    V1=np.zeros(Nx+1)
    x_cor = np.linspace(lower_x, lx, Nx + 1)
    ctr=0 #marker for the value in step_arr
    #if alpha<=0.5:print("Stability can be a problem")
    #SET fig size
    plt.figure(figsize=(10,10))
    for i in range(0,Nx+1):
        if lower_x + (hx * i) == 1:#1 as inital 300C at length 1
            V0[i]=300
        else:
            V0[i]=0
        x_cor[i]=(lower_x + hx * i)
    plt.plot(x_cor, V0, label=0)
    #Matrix mult for sparse when only some are multiplied
    for j in range(0,tot_steps+1):#1000 is number of steps taken
        for i in range(0,Nx+1):

            if i==0:
                V1[i]=(1-2*alpha)*V0[i]+alpha*V0[i+1]
            elif i==Nx:
                V1[i]=(1-2*alpha)*V0[i]+alpha*V0[i-1]
            else:
                V1[i]=(1-2*alpha)*V0[i]+alpha*V0[i-1]+alpha*V0[i+1]
        for k in range(0,Nx+1):#Equating array V0 to V1
            V0[k]=V1[k]
        if j==step_arr[ctr]:
            plt.plot(x_cor,V1,label=step_arr[ctr])
            #print(V0[50])
            ctr=ctr+1
    plt.legend()
    return None

def Eigen_solve(A, k_max, tol, guess,r):
    #tol = 10^(-r)
    na = len(A)
    b_k=[]
    lambda1=0
    if str(guess)=='random':
        for i in range(0,na):
            R=[]
            r=toolsar.rnum(i)
            R.append(r)
            b_k.append(R)
    else:
        b_k = guess
    #print('\n Guess vector:',b_k)
    
    e_new=1
    e_old=0
    i=0
    while (abs(e_new - e_old) >= tol) and (i in range(0,k_max)):
        e_old=lambda1
        b_k = toolsar.crossmat(A,b_k)
        b_k_n = 0
        for j in range(0,len(b_k)):
            b_k_n += (b_k[j][0])**2
        
        for j in range(0,len(b_k)):
            b_k[j][0] = b_k[j][0] / math.sqrt(b_k_n)
        
        i=i+1
        
        #finding the value of b_k*Ab_k
        x1=toolsar.matrixtranspose(b_k)
        x2=toolsar.crossmat(x1,A)
        x3=toolsar.crossmat(x2,b_k)
        #finding value of b_k*b_k
        v=toolsar.matrixtranspose(b_k)
        v1=toolsar.crossmat(v,b_k)
        #lambda1 = b_k*Ab_k / b_k*b_k
        lambda1=x3[0][0]/v1[0][0]
        
        e_new=lambda1
    
    for l in range(len(A)):
        A[l] = [round(ele,r) for ele in A[l]]
        b_k[l][0] = round(b_k[l][0],r)
    print('\n Matrix A:',A,'\n \n Eignvector:',b_k,f'\n \n Eigenvalue:{lambda1:1.{r}f}','\n \n', i,'iterations.')
    return b_k,lambda1

#--------------------------------------------------------#
#-------------Fixed Point Iteration Method---------------#
#--------------------------------------------------------#
def fixed_point_iteration(x0 , g , tol):
    """This function is used to find the root of a given function using fixed point iteration method.

    Args:
        x0 (int or float): initial guess of the root.
        g (function): Function g(x) = x.
        tol (float): tolerance.

    Returns:
        float: root of the given function to the specified tolerance.
    """
    x = x0; iter = 0
    while True:
        x = g(x)
        iter += 1
        if abs(x - x0) < tol:
            break
        x0 = x
    #rounding off the value of x to given tolerance
    x = round(x, int(abs(math.log10(tol))))
    return x, iter
#--------------------------------------------------------#

#--------------------------------------------------------#
#--------------------Simpson's Rule----------------------#
#--------------------------------------------------------#
def simpson_int(fn,x0,x2,n):
    """This function is used to find the integration of a given function using Simpson's rule.

    Args:
        fn (function): function to be integrated.
        x0 (float): lower limit of the integration.
        x2 (float): upper limit of the integration.
        n (int): number of intervals.

    Returns:
        float : integration of the given function.
    """
    h = abs(x2-x0)/n
    add = 0
    xi = x0
    for i in range(n+1):
        if i == 0 or i == n:
            add += (h/3)*fn(xi)
        elif i%2 == 0:
            add += (2*h/3)*fn(xi)
        else:
            add += (4*h/3)*fn(xi)
        xi += h
    return add
#--------------------------------------------------------#

#--------------------------------------------------------#
#-------------------Gaussian Quadrature------------------#
#--------------------------------------------------------#
def gauss_nodes_and_weights(n, a, b):
    """This function is used to find the nodes and weights for Gaussian Quadrature.

    Args:
        n (int): number of intervals.
        a (float): lower limit of the integration.
        b (float): upper limit of the integration.

    Returns:
        float : nodes and weights for Gaussian Quadrature.
    """
    x = []
    w = []
    if n == 1:
        x.append(0.0)
        w.append(2.0)
        return x, w
    for i in range(n):
        x0 = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        x1 = 0
        while abs(x0 - x1) > 1e-10:
            p0 = 1.0
            p1 = 0.0
            for j in range(n):
                p2 = p1
                p1 = p0
                p0 = ((2 * j + 1) * x0 * p1 - j * p2) / (j + 1)
            pp = n * (x0 * p0 - p1) / (x0 * x0 - 1)
            x1 = x0
            x0 = x1 - p0 / pp
        x.append(x0)
        w.append(1 / ((1 - x0 * x0) * pp * pp))
    return x, w

def gauss_quadrature(fn, a, b, n):
    """This function is used to find the integration of a given function using Gaussian Quadrature.

    Args:
        fn (function): function to be integrated.
        a (float): lower limit of the integration.
        b (float): upper limit of the integration.
        n (int): number of intervals.

    Returns:
        float : integration of the given function.
    """
    x, w = gauss_nodes_and_weights(n, a, b)
    integral = 0.0
    for i in range(n):
        integral += w[i] * fn(x[i])
    return integral
#--------------------------------------------------------#

#--------------------------------------------------------#
#-------------------------Rk4----------------------------#
#--------------------------------------------------------#
def rk4_step(f, x, y, h):
    """Performs a single step of the fourth-order Runge-Kutta method.

    Args:
        f (function): The derivative y' in terms of x and y.
        x (float): The current x value.
        y (float): The current y value.
        h (float): The step size.

    Returns:
        float: The updated y value after the RK4 step.
    """
    k1 = h * f(x, y)
    k2 = h * f(x + h/2, y + k1/2)
    k3 = h * f(x + h/2, y + k2/2)
    k4 = h * f(x + h, y + k3)
    return y + (k1 + 2*k2 + 2*k3 + k4) / 6

def solve_ode_rk4(f, initial_x, initial_y, interval_size, num_steps):
    """Solves a first-order ordinary differential equation using the fourth-order Runge-Kutta method.
    
    Parameters:
    - f: The derivative y' in terms of x and y.
    - initial_x: The initial value of the independent variable.
    - initial_y: The initial value of the dependent variable.
    - interval_size: The size of the interval between each step.
    - num_steps: The number of steps to take.
    
    Returns:
    - solution: A list of tuples representing the solution, where each tuple contains the x and y values.
    """
    solution = [(initial_x, initial_y)]
    x = initial_x
    y = initial_y
    for _ in range(num_steps):
        y = rk4_step(f, x, y, interval_size)
        x += interval_size
        solution.append((x, y))
    return solution

#--------------------------------------------------------#

#--------------------------------------------------------#
#-------------------Crank-Nicolson Method----------------#
#--------------------------------------------------------#
def crank_nicolson(g, L, n, T, dt):
    """Solve the 1D heat equation using the Crank-Nicolson method.

    Parameters:
    g (function): The initial condition function g(x).
    L (float): The length of the domain.
    n (int): The number of grid points.
    T (float): The total time.
    dt (float): The time step size.

    Returns:
    tuple: A tuple containing the solution at each time step and the grid points.
    """

    dx = L / (n + 1)
    alpha = dt / (dx**2)
    A = np.zeros((n, n))
    B = np.zeros((n, n))

    # Create tridiagonal matrices A and B
    for i in range(n):
        A[i, i] = 2 + 2 * alpha
        B[i, i] = 2 - 2 * alpha
        for j in range(n):
            if j == i + 1 or j == i - 1:
                A[i, j] = -alpha
                B[i, j] = alpha

    x = [0]
    for i in range(n - 1):
        x.append(x[i] + dx)

    # Initialize initial values of the vector v0 using g(x)
    v0 = np.array([g(xi) for xi in x])
    v0[-1] = v0[0]

    v = v0.copy()
    solution_at_each_time = [v0.copy()]

    # Time-stepping using Crank-Nicolson method
    for _ in range(int(T / dt)):
        C = np.matmul(np.linalg.inv(A), B)
        v = np.matmul(C, v)
        solution_at_each_time.append(v.copy())

    return solution_at_each_time, x
#--------------------------------------------------------#

def QR_eigenfind(A,tolerance = 1e-6):
    A = np.array(A)
    copy_A = np.copy(A)
    Q,R = QR_factorize(A)
    A = np.matmul(R,Q)
    i=1
    while np.linalg.norm(A-copy_A)>tolerance:
        copy_A = np.copy(A)
        Q,R = QR_factorize(A)
        A = np.matmul(R,Q)
        i+=1
    return np.diag(A),i

def QR_factorize(A):
    A = np.array(A) if type(A) != np.ndarray else A
    Q = np.zeros(A.shape)
    R = np.zeros(A.shape)
    for i in range(A.shape[1]):
        u_i = A[:,i]
        sum = 0
        for j in range(i):
            sum += np.dot(A[:,i],Q[:,j])*Q[:,j]
        u_i = u_i - sum
        Q[:,i] = u_i/np.linalg.norm(u_i)
        for j in range(i+1):
            R[j,i] = np.dot(A[:,i],Q[:,j])
            
    return Q,R

def polynomial_fit(xlist: list,ylist: list,sigma_list: list,degree: int,tol=1e-6):
    xlist = np.array(xlist)
    ylist = np.array(ylist)
    sigma_list = np.array(sigma_list)
    A_matrix = np.zeros((degree+1,degree+1))

    for i in range(degree+1):
        for j in range(degree+1):
            A_matrix[i][j] = np.sum((xlist**(i+j))/(sigma_list**2))
    B_matrix = np.zeros(degree+1)
    for i in range(degree+1):
        B_matrix[i] = np.sum((ylist*(xlist**i))/(sigma_list**2))
    # a = Gauss_seidel_solve(A_matrix.tolist(),B_matrix.tolist(),T=tol)
    a = np.linalg.solve(A_matrix,B_matrix)    
    return a,A_matrix

def poly_fn(x,coefflist):
    sum = 0
    for i in range(len(coefflist)):
        sum += coefflist[i]*x**i
    return sum

def polynomial_fit_mod_chebyshev(xlist: list,ylist: list,sigma_list: list,degree: int):
    # Defining the modified chebyshev polynomial
    def modified_chebyshev_polynomial(x,degree):
        def chebyshev_polynomial(x,degree):
            if degree == 0:
                return 1
            elif degree == 1:
                return x
            else:
                return 2*x*chebyshev_polynomial(x,degree-1) - chebyshev_polynomial(x,degree-2)
        return chebyshev_polynomial(2*x - 1,degree)
    xlist = np.array(xlist)
    ylist = np.array(ylist)
    sigma_list = np.array(sigma_list)
    A_matrix = np.zeros((degree+1,degree+1))

    for i in range(degree+1):
        for j in range(degree+1):
            # Replace the polynomial with the modified chebyshev polynomial
            A_matrix[i][j] = np.sum((modified_chebyshev_polynomial(xlist,i)*modified_chebyshev_polynomial(xlist,j))/(sigma_list**2))
    B_matrix = np.zeros(degree+1)
    for i in range(degree+1):
        B_matrix[i] = np.sum((ylist*(modified_chebyshev_polynomial(xlist,i)))/(sigma_list**2))
    a = np.linalg.solve(A_matrix,B_matrix)    
    return a,A_matrix

def modified_chebyshev_polynomial(x,degree):
    def chebyshev_polynomial(x,degree):
        if degree == 0:
            return 1
        elif degree == 1:
            return x
        else:
            return 2*x*chebyshev_polynomial(x,degree-1) - chebyshev_polynomial(x,degree-2)
    return chebyshev_polynomial(2*x - 1,degree)


def poly_fn_mod(x,coefflist):
    sum = 0
    for i in range(len(coefflist)):
        sum += coefflist[i]*modified_chebyshev_polynomial(x,i)
    return sum   

#####################################################################################
#                            Random Number Generation                             
#####################################################################################
class rng():
    def __init__(self,seed, a = 1103515245, c = 12345 ,m = 32768):
        # initiation of data input
        self.term = seed
        self.a = a
        self.c = c
        self.m = m
    def gen(self):
        # generates a random number
        self.term = (((self.a * self.term) + self.c) % self.m)
        return self.term/self.m
    def genlist(self,length):
        # returns a list of 'n' random numbers in the range (0,1) where 'n' is 'length'.
        RNs = []
        for i in range(length):
            self.term = (((self.a * self.term) + self.c) % self.m)
            RNs.append(self.term / self.m)
        return RNs
    
def monte_carlo_integrate(f: float,a: float,b: float,N: int,seed: int,multiplier=1103515245,m=32768,c=12345):
    '''
    # Monte Carlo Integration
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - N: Number of random numbers to be generated
    - seed: Seed for the random number generator
    ## Returns
    - F: The value of the integral
    '''
    p=rng(seed,m=m,c=c,a=multiplier)
    F=0
    for i in range(N):
        k=p.gen()
        k=((b-a)*(k))+a
        F+=((b-a)*f(k))/N   
    return F   

def monte_carlo_error(f: float,a: float,b: float,N: int,seed: int):
    '''
    # Monte Carlo Integration
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - N: Number of random numbers to be generated
    - seed: Seed for the random number generator
    '''
    rn=rng(seed)
    F=0
    F1=0
    for i in range(N):
        p=rn.gen()
        p=((b-a)*(p/32768))+a
        F+=f(p)
        F1+=pow(f(p),2)  
    return (F1/N)-pow((F/N),2) 
