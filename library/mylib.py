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

# function to print a matrix aesthetically
def Print_Matrix(a):[print((' \t').join(str(k) for k in a[i]),'\n') for i in range(len(a))]
    
#Function to read the arrays from an external .txt file
def read_arr(path):
    contents = {}
    with open(path,'r') as f:
        data = f.readlines()
    #print(data)
    for line in data:
        lin_lst = line.strip().split()
        #print(lin_lst)
        if lin_lst[0] == '#':
            name = lin_lst[1]
            contents[name] = []
        else:
            if 'matrix' in name:
                contents[name].append([float(i) for i in lin_lst])
            elif 'vector' in name:
                contents[name] = [float(i) for i in lin_lst]
    return contents

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
def LU_solve(A,B,r=6):
    L,U = LU_decompose(A)
    y = forward_sub(L,B)
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
        print("You have not assumed right a and b\n")
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
        return zl, yl, xl
    else:
        return z_i, y_i

def boundary_RK4_dyx(y,z,x):
    f = z
    return f

def boundary_RK4_dzx(y,z,x):
    f = 0.01*(y-20)
    return f

import numpy as np
import toolsar, math

def PDE_Solve(lx,Nx,lt,Nt,lower_x,tot_steps):
    step_arr = [10, 20, 50, 100, 200, 500, 1000, tot_steps]
    hx=(lx/Nx)
    ht=(lt/Nt)
    alpha=ht/(hx)**2
    V0=np.zeros(Nx+1)
    V1=np.zeros(Nx+1)
    x_cor = np.linspace(lower_x, lx, Nx + 1)
    ctr=0 #marker for the value in step_arr
    #if alpha<=0.5:print("Stability can be a problem")
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
