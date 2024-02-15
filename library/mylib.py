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
def Simpson_int(fn,x0,x2,n):
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
def rk4_step(x, y, h):
    """Performs a single step of the fourth-order Runge-Kutta method.

    Args:
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

def solve_ode_rk4(initial_x, initial_y, interval_size, num_steps):
    """ Solves a first-order ordinary differential equation using the fourth-order Runge-Kutta method.
    
    Parameters:
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
        y = rk4_step(x, y, interval_size)
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
