#Random Number Generator
class pRNG():
    '''My class to store seed value and different types of pRNG'''
    seed = 10
    def LCG(low=0,high=1):
        '''Gives a random number between low and high'''
        a = 1103515245
        c = 12345
        m = 32768
        x1 = ((a*(pRNG.seed)+c)%m)/m
        pRNG.seed = x1 #changing the seed value every time the rng is called
        val = low + abs(high-low)*x1
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
def GJ_Inverse(a,r_factor = 5):
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

#Functions for root finding

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
        print(f"Iteration no: {n} \troot -> {x:{1}.{e}f}")
    print(f"The root converges to {x:{1}.{e}f} and the value of function is {fn(x):{1}.{e}f}")
    print("\n Total no of iterations = ",n)
    #return x


#Laguerre's Method

# Function to find the value of a polynomial at a given 'x'
# using its list of coefficients
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

#Data Interpolation and Fitting

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
