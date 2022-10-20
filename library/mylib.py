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
def Gauss_Jordan_Solve(a,b):
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
  A = [round(ele[-1],3) for ele in A]
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
def determinant(A):
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
  Print_Matrix(a)
  det = (-1)**count
  for l in range(n):
    det *= a[l][l]
  return round(det,4)

#L-U Decomposition

#Function to decompose a square matrix into lower and upper matrix
def LU_decompose(A):
  n = len(A)
  U = [A[i][:] for i in range(n)]
  L = Identity(n)
  for i in range(n):
    if swap(U,i)==None:
      print("Solution is not possible!")
      return None
    for j in range(i+1,n):
      L[j][i] = U[j][i]/U[i][i]
      U[j] = [U[j][k]-L[j][i]*U[i][k] for k in range(n)]
    L = [[round(x,6) for x in ele] for ele in L]
    U = [[round(x,6) for x in ele] for ele in U]
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
    x = [round(i,3) for i in x]
  return x

#Function to solve a system of linear equation using LU-Decomposition
def LU_solve(A,B):
  L,U = LU_decompose(A)
  y = forward_sub(L,B)
  x = backward_sub(U,y)
  return x

#Cholesky Decomposition

#Function to check whether the given matrix is symmetric or not
def isSymmetric(a):
  n = len(a)
  try:
    for i in range(n):
      for j in range(n):
        if a[i][j] != a[j][i]:
          return False
    return True
  except:
    print("matrix is not a square matrix")
    
#Function to take a transpose of given matrix
def transpose(a):
  if type(a[0]) == int or type(a[0]) == float or type(a[0]) == complex: 
    return [[i] for i in a]
  elif type(a[0]) == list and len(a[0]) == 1:
    return [a[i][0] for i in range(len(a))]
  tp = [len(a)*[0] for _ in range(len(a[0]))]
  for i in range(len(a)):
    for j in range(len(a[0])):
      tp[i][j] = a[j][i]
  return tp

#Cholesky decomposition function returning an Upper matrix
def Cholesky_decompose(a):
  if isSymmetric(a) == False:
    print('Matrix is not symmetric')
    return None
  final = [len(a)*[0] for _ in range(len(a))]
  for i in range(len(a)):
    for j in range(i+1):
      add = sum(final[i][k]*final[j][k] for k in range(j))
      if i == j:
        final[i][j] = (a[i][i]-add)**0.5
      else:
        final[i][j] = 1/final[j][j]*(a[i][j]-add)
  #final = [[round(x,3) for x in ele] for ele in final]
  return final

#Function to solve system f equations using cholesky decomposition
def Cholesky_solve(A,B):
  U = Cholesky_decompose(A)
  L = transpose(U)
  for i in range(len(L)):
    L[i][i] = 1
  y = forward_sub(L,B)
  x = backward_sub(U,y)
  return x

#Jacobi Method

#Function to check whether the given matrix is diagonally dominant
def isDiaDominant(A):
  for i in range(len(A)):
    if not abs(A[i][i]) >= sum(abs(A[i][j]) for j in range(len(A)) if j!=i):
      return False
  return True

# Swap two given rows of a matrix
def swap1(a,r1,r2):
    for j in range(len(a)):
        a[r1][j],a[r2][j] = a[r2][j],a[r1][j]
    return a

# make a matrix diagonally dominant 
def diagDM (a,b):
    n = len(a)
    if isDiaDominant(a) != True:
        sum = 0 
        for i in range(n):
            for j in range(n):
                if j != i:
                    sum += abs(a[i][j])
            if abs(a[i][i])<sum:
                for k in range(i,n):
                    a = swap1(a,i,k)
                    b[i],b[k] = b[k],b[i]
                    s = 0 
                    for m in range(n):
                        if m != i:
                            s += abs(a[i][m])
                    if a[i][i] > s:
                        break
                    else:
                        continue
            else:
                sum = 0
                continue
    else:
        return a,b

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
      condition = False not in [abs(xi[k]-xj[k])>10**(-1*e) for k in range(n)]
      #print(f"Iteration number: {count}, {[f'{i:{1}.{e}}' for i in xi]}\n")
    xi = [round(i,e) for i in xi]
    print(f"No of Iterations: {count}\n")
    return xi

#Functions for root finding

#Function to find root using Bisection method
def Bisection(fn,a,b,e=6):
    if (fn(a) * fn(b) >= 0):
        print("You have not assumed right a and b\n")
        return None
    c = a
    n = 0
    while (b-a) >= 10**-e:
        n += 1
        c = (a+b)/2
        if abs(fn(c)) <= 10**-e:
            break
        if (fn(c)*fn(a)<0):
            b = c
        else:
            a = c
    print(f"The root in the given interval is {c:{1}.{e}} and the value of function is {fn(c):{1}.{e}f}")
    print("No. of iterations = ",n)
    return round(c,e)

#Function to show bracketing
def Bracketing(fn,a,b,e=4):
    if (fn(a) * fn(b) >= 0):
        print("You have not assumed right a and b\n")
        return None
    c = a
    n = 0
    while (b-a) >= 10**-e:
        n+=1
        c = (a+b)/2
        if abs(fn(c)) <= 10**-e:
            break
        if (fn(c)*fn(a)<0):
            b = c
            print(f"N: {n} \t\tBracketing: \t({a:{1}.{e}},{c:{1}.{e}})")
        else:
            a = c
            print(f"N: {n} \t\tBracketing: \t({c:{1}.{e}},{b:{1}.{e}})")
    print("Total no of iterations = ",n-1)
    
#Function to find root using Regula-Falsi method
def Regula_falsi(fn,a,b,e=6):
    if (fn(a) * fn(b) >= 0):
        print("You have not assumed right a and b\n")
        return None
    n = 0
    c=a
    c1 = c-1
    while abs(c-c1) > 10**-e or abs(fn(c)) > 10**-e:
        n += 1
        c = (a * fn(b) - b * fn(a))/ (fn(b) - fn(a))
        if abs(fn(c)) <= 10**-e:
            break
        if (fn(c)*fn(a)<0):
            b = c
        elif (fn(c)*fn(b)<0):
            a = c
        c1 = c
    print(f"The root in the given interval is {c:{1}.{e}} and the value of function is {fn(c):{1}.{e}f}")
    print("Total no of iterations = ",n)
    return round(c,e)

#Function to find root using Newton_raphson method
def Newton_raphson(fn,d_fn,x = 0.5,e=6):
    '''fn: the function of which we want to find the root,
       d_fn: the derivative of the function
       x: initial guess for the root'''
    h = fn(x)/d_fn(x)
    n = 0
    while abs(h)>10**-e:
        n+=1
        h = fn(x)/d_fn(x)
        x = x-h
    #print(f"The root is {x:{1}.{e}f} and the value of function is {fn(x):{1}.{e}f}")
    print("Total no of iterations = ",n)
    return x

