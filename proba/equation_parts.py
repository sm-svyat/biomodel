from math import exp, log10

def psih(u, v, r, s, p, rho):
    '''
    psih(u, v, r, p, rho)
    nonlinear part in 0th and 8th equations
    '''
    c20h = p['cc20h']
    c23h = p['cc23h']
    cH01 = p['ccH01']
    c21h = p['cc21h']
    c22h = p['cc22h']
    c24h = p['cc24h']
    cc5r  = p['cc5r']
    cc6r  = p['cc6r']
    r0    = p['r0']

    return (v-c20h*rho*u + c23h*cH01*u*v)/(c21h + v + c22h*rho*u + c24h*cH01*u*v)*(cc5r + cc6r*r0*r)


def psig(u, p):
    '''
    psig(u, p)
    nonlinear part in 0th and 8th equations
    '''
    c4s = p['cc4s']
    c40s = p['cc40s']
    return  ((u/(c4s*c40s)))


def psicc(m, s, u, o, r, p):
    '''
    psicc(m, s, u, o, r, p)
    nonlinear part in equations 
    '''
    c1s   = p['cc1s']
    c2s   = p['cc2s']
    c3s   = p['cc3s']
    c6s   = p['cc6s']
    c40s  = p['cc40s']
    c41s  = p['cc41s']
    c42s  = p['cc42s']
    csmin = p['csmin']
    csmax = p['csmax']
    cc7r  = p['cc7r']
    cc8r  = p['cc8r']
    r0    = p['r0']

    return (m*s)/(m + c2s)/(s + c3s)*((csmin - csmax)/(1 + exp((c1s - u)/c6s)) + csmax)*((c40s)/(c40s + c41s*(1 + o/c42s)))*(cc7r + cc8r*r0*r)


def psi2(u, v, rho, p):
    '''
    nonlinear part in 7th and 8th equations 
    '''
    c12 = p['cc12']
    c13 = p['cc13']

    return (v-rho*u)/(c12 + c13*(v+rho*u))


def psis(u, p):
    '''
    psis(u, p)
    nonlinear part in 0th and 8th equations
    '''
    u0  = p['u0']
    c7s = p['cc7s']
    c8s = p['cc8s']

    return (u - u0)/(c7s + c8s*(u+u0))


def y2_exp(x, **params):
    '''
    y2_exp(x, **params)
    doing samething I have no idea about
    '''
    a  = params['a']
    b  = params['b']
    k  = params['k']
    x0 = params['x0']
    C  = (b - a) / 2
    D  = a
    B  = (b+a) / 2


    xsh = x - x0

    if xsh >= 0:
        y = C*(2 - exp(-k*xsh)) + D
        y_t = C*k*xsh + B
    else:
        y = C*exp(k*xsh) + D
        y_t = C*k*xsh + B
    return y


def Light2(v, p):
    '''
    Light2(v, p)
    nonlinear illumination in 3th and 6th equations 
    '''    
    Light2_pHmin   = p['Light2_pHmin']
    Light2_pHmax   = p['Light2_pHmax']
    Light2_pH0     = p['Light2_pH0']
    Light2_Lmin    = p['Light2_Lmin']
    Light2_Lmax    = p['Light2_Lmax']
    Light2_tangent = p['Light2_tangent']

    x = 6 - log10(v)

    return y2_exp(x, a=Light2_Lmin, b=Light2_Lmax,k=Light2_tangent, x0=Light2_pH0,xmin=Light2_pHmin,xmax=Light2_pHmax)


def tau(y, z, v, p):
    '''
    tau(y, z, v, p)
    common nonlinear part in f2 and g2
    but I have no idea what do f2 and g2 mean
    '''
    k1 = p['kk1']
    k2 = p['kk2']
    k3 = p['kk3']
    return (1/(1-z))/k1+ 1/y/k2 + (1/k3)*(1+v)


def alpha(u, p):
    '''
    buffer function pHout
    :param u: переменная. конецентрация протонов в строме
    :param p: словарь с константами
    :return:
    '''
    return 1/p['aa1']
    #return 1 + (p['aa1'])/(p['aa2'] + p['aa3']*u)/(p['aa2'] + p['aa3']*u)


def beta(v, p):
    '''
    buffer function pHin
    :param v:
    :param p:
    :return:
    '''
    return 1/p['bb1']
    #return 1 + (p['bb1'])/(p['bb2'] + p['bb3']*v)/(p['bb2'] + p['bb3']*v)