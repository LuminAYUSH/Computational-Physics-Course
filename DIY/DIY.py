def adaptive_step_rk(h, t, t_lim, w, func, epsilon):
    """
    h : initial step size
    t : time parameter
    w : initial value for solution curve y
    func : a function of t and y
    epsilon : tolerance value used for the solution
    
    T : list of values of t
    Y : list of values computed for y for corresponding values of t
    """
    T = []
    Y = []
    
    while t < t_lim:
        h = float(min(h, t_lim-t))
        
        #increments for the RK algorithm
        k1 = h*func(t, w)
        k2 = h*func(t + (1/4)*h, w + (1/4)*k1)
        k3 = h*func(t + (3/8)*h, w + (3/32)*k1 + (9/32)*k2)
        k4 = h*func(t + (12/13)*h, w + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
        k5 = h*func(t + h, w + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        k6 = h*func(t + (1/2)*h, w - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
        
        #summations of weights for increments
        w1 = w + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5
        w2 = w + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
        
        if w1-w2 == 0:
            TE = 0.00000001   #TE set to minimum possible float value
        else:
            TE = abs(w1-w2)/h #truncation error calculated
            
        delta = 0.95*(epsilon/TE)**(1/4)
        
        if TE <= epsilon: #step size revised while moving to next step
            t = t+h
            w = w1
            T.append(t)
            Y.append(w)
            h = delta*h
        else:            #step size revised & calculation is repeated
            h = delta*h
            
    return T, Y

def adaptive_step_rk_2var(h, t, t_lim, w, l, func1, func2, epsilon):
    """
    h : initial step size
    t : initial value for time parameter
    t_lim : final value for time parameter
    w : initial value for solution curve x
    l : initial value for solution curve y
    func1 : a function of t, x, and y with solution curve x
    func2 : a function of t, x, and y with solution curve y
    epsilon : tolerance value used for the solution
    
    T : list of values of t
    Y : list of values computed for y for corresponding values of t
    """
    T = []
    X = []
    Y = []
    
    while t < t_lim:
        h = float(min(h, t_lim-t)) #smallest possible step size is chosen near the endpoint
        
        #increments for the RK algorithm for func1
        k1_1 = h*func1(t, w, l)
        k1_2 = h*func1(t + (1/4)*h, w + (1/4)*k1_1, l)
        k1_3 = h*func1(t + (3/8)*h, w + (3/32)*k1_1 + (9/32)*k1_2, l)
        k1_4 = h*func1(t + (12/13)*h, w + (1932/2197)*k1_1 - (7200/2197)*k1_2 + (7296/2197)*k1_3, l)
        k1_5 = h*func1(t + h, w + (439/216)*k1_1 - 8*k1_2 + (3680/513)*k1_3 - (845/4104)*k1_4, l)
        k1_6 = h*func1(t + (1/2)*h, w - (8/27)*k1_1 + 2*k1_2 - (3544/2565)*k1_3 + (1859/4104)*k1_4 - (11/40)*k1_5, l)
        
        #increments for the RK algorithm for func2
        k2_1 = h*func2(t, w, l)
        k2_2 = h*func2(t + (1/4)*h, w, l + (1/4)*k2_1)
        k2_3 = h*func2(t + (3/8)*h, w, l + (3/32)*k2_1 + (9/32)*k2_2)
        k2_4 = h*func2(t + (12/13)*h, w, l + (1932/2197)*k2_1 - (7200/2197)*k2_2 + (7296/2197)*k2_3)
        k2_5 = h*func2(t + h, w, l + (439/216)*k2_1 - 8*k2_2 + (3680/513)*k2_3 - (845/4104)*k2_4)
        k2_6 = h*func2(t + (1/2)*h, w, l - (8/27)*k2_1 + 2*k2_2 - (3544/2565)*k2_3 + (1859/4104)*k2_4 - (11/40)*k2_5)
        
        #summations of weights for increments for func1
        w1 = w + (25/216)*k1_1 + (1408/2565)*k1_3 + (2197/4104)*k1_4 - (1/5)*k1_5
        w2 = w + (16/135)*k1_1 + (6656/12825)*k1_3 + (28561/56430)*k1_4 - (9/50)*k1_5 + (2/55)*k1_6
        
        #summations of weights for increments for func2
        l1 = l + (25/216)*k2_1 + (1408/2565)*k2_3 + (2197/4104)*k2_4 - (1/5)*k2_5
        l2 = l + (16/135)*k2_1 + (6656/12825)*k2_3 + (28561/56430)*k2_4 - (9/50)*k2_5 + (2/55)*k2_6
        
        if w1-w2 == 0:
            TE_1 = 0.00000001   #TE set to minimum possible float value
        else:
            TE_1 = abs(w1-w2)/h #truncation error calculated
            
        if l1-l2 == 0:
            TE_2 = 0.00000001   #TE set to minimum possible float value
        else:
            TE_2 = abs(l1-l2)/h #truncation error calculated
            
        delta1 = 0.9*(epsilon/TE_1)**(1/4)
        delta2 = 0.9*(epsilon/TE_2)**(1/4)
        
        if delta1 < delta2: #the smaller delta chosen to ensure maximum smoothness
            delta = delta1
        else:
            delta = delta2
        
        if TE_1 <= epsilon and TE_2 <= epsilon: #step size revised while moving to next step
            t = t+h
            w = w1
            l = l1
            T.append(t)
            X.append(w)
            Y.append(l)
            h = delta*h
        else:            #step size revised and calculation is repeated
            h = delta*h
            
    return T, X, Y
