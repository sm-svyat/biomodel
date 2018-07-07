#!/usr/bin/python3
import sys
sys.path.append('~/git/pyproba/server/proba')
from equation_parts import psih, psig, psicc, psi2, psis, Light2, tau, alpha, beta
from math import log10, exp


def read_constants(CONST_PATH):
    '''
    Функция читает константы из ini-файла матлабовской модели
    :param CONST_PATH:
    :return:
    '''
    constants={}
    with open(CONST_PATH, 'r') as const_file:
        for row in const_file:
            try:
                const_name, const_values = row.split('#')
            except ValueError:
                continue
            const_value = const_values.split('%')
            constants[const_name.split(' ')[0]]=float(const_value[0])
    return constants


def equcoeff(p):
    '''
    Функция переопеределения параметров
    :param p:
     equ_switch = [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]:
                  [s, w, o, q, x, y, z, u, v, n, m, r]

    equZ_switch = [1,    1,    1,    1,    1,    1,    1,    1,    1]
                  [x0,   x1,  1-a0, 1-b0, 1-a1, 1-b1, 1-fx, 1-f,  mv]
    :return:
    '''

    #Константы обратных реакций для детализированной ФС1
    p['ta0'] = p['ca0'] * exp(-p['dE_a0'] / 26) # kT = 26mv
    p['ta1'] = p['ca1'] * exp(-p['dE_a1'] / 26)
    p['tfa'] = p['cfa'] * exp(-p['dE_fa'] / 26)
    p['tff'] = p['cff'] * exp(-p['dE_ff'] / 26)
    p['tmv'] = p['cmv'] * exp(-p['dE_mv'] / 26)
    p['tb0'] = p['cb0'] * exp(-p['dE_b0'] / 26)
    p['tb1'] = p['cb1'] * exp(-p['dE_b1'] / 26)
    p['tfb'] = p['cfb'] * exp(-p['dE_fb'] / 26)
    p['txx'] = p['cxx'] * exp(-p['dE_xx'] / 26)
    #p['txx'] = 1

    # artifisial
    p['tx1'] = 0.5
    p['tx2'] = 0.5
    p['tx3'] = 0.5
    p['tx4'] = 0.05
    p['tx5'] = 0.01 #0.01
    p['tx6'] = 0.05 #0.05

    # Определение списка подключенных уравнений equ_switch (стандартная модель), equZ_switch (детализация ФС1)
    #p['equ_switch'] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    p['equZ_switch'] = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    p['equ_switch'] = [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    #p['equZ_switch'] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    return p


def boundary_conditions(y):
    '''
    The function takes into account boundary conditions of the model
    :param y:
    :return y:
    '''
    u_bnd = 1 - 1e-7
    l_bnd = 1e-7
    for k in range(len(y)):
        if k != 7 and k != 8:
            if y[k] > u_bnd:
                y[k] = u_bnd
            elif y[k] < l_bnd:
                y[k] = l_bnd
        else:                   #следует проверить целесобразность этого условия #TODO
            if y[k] < l_bnd:
                y[k] = l_bnd
    return y


def equinit(y, p):
    '''
    setting equations  coefficients constants(no arguments)
    :param y:
    :param equ_switch = [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]:
                        [s, w, o, q, x, y, z, u, v, n, m, r]

          equZ_switch = [1,    1,    1,    1,    1,    1,    1,    1,    1]
                        [x0,   x1,  1-a0, 1-b0, 1-a1, 1-b1, 1-fx, 1-f,  mv]
    :return mass(y):
    '''
    try:
        equ_switch = p['equ_switch']
        equZ_switch = p['equZ_switch']
    except KeyError:
        print('Do not specify a list of enabled equations')
        p['equ_switch'] = [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        p['equZ_switch'] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        equ_switch = p['equ_switch']
        equZ_switch = p['equZ_switch']
    y = [i[0] * i[1] for i in zip(equ_switch+equZ_switch, y)]
    return mass(y, p)


def mass(y, p):
    # Mass matrix function
    m = [1, 1, 1, 1, 1, 1, 1, alpha(y[7], p), beta(y[8], p), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    return [i[0]*i[1] for i in zip(m, y)]


def f(t, vars, p):
    '''
    The equations' system of the model
    :param t:
    :param vars:
    :return equinit(wrs):
    '''
    #сдвинул индексы на 1 для облегчения читаемости системы уравнений
    y = [t, vars[0], vars[1], vars[2], vars[3], vars[4], vars[5], vars[6], vars[7], vars[8], vars[9], vars[10], vars[11]]
    z = [t] + [vars[i] for i in range(12, len(vars)) ]

    rho = (10**(-p['ccH01'])/(10**(-p['ccH1'])))
    L1 = p['cL1']
    L2 = p['cL2']

    # 1    2    3    4    5    6     7    8    9     10     11     12
    # s    w    o    q    x    y     z    u    v     n      m      r

    wr1 = p['cc2h'] * p['cc11h'] * (1 - y[1]) * psih(y[8], y[9], y[12], y[1], p=p, rho=rho) - p['cc1h'] * p['cc3h'] * p['n0'] * (1 + p['cca'] * psig(y[3], p=p)) * psicc(y[11], y[1], y[8], y[3], y[12], p=p) - p['cc4h'] * y[1]

    #wr2 = (p['cc8']*p['z0']*y[8]*y[7] + p['cc34']*p['n0']*y[10] + p['cc1r']*y[8]*p['r0']*(1 - y[12]) + p['cc6']*p['o0']*y[3])*(1 - y[2]) - p['cL1']*p['cc7']*(1 - y[5])*y[2]
    wr2 = (p['cc8'] * p['z0'] * y[8] * y[7] + p['cc34'] * p['n0'] * y[10] + p['cc1r'] * y[8] * p['r0'] * (1 - y[12]) + p['cc6'] * p['o0'] * y[3]) * (1 - y[2]) - p['cL1'] * p['cc7'] * (z[8]) * y[2] # Detailed PSI

    wr3 = 0.25 * p['cc30'] * p['q0'] / p['o0'] * y[4] - 0.5 * y[3] * (p['cc31'] * p['z0'] * (1 - y[7]) + p['cc32'] * p['y_0'] * (1 - y[6]) + p['cc6'] * p['w0'] * (1 - y[2]) + 0 * 1 * y[12]) - 0.0000 * (y[3] - 0.5)

    wr4 = p['cL2'] * Light2(y[9], p) * p['cc10'] * p['z0'] * y[8] * y[7] * (1 - y[4]) - p['cc30'] * y[4]

    wr5 = p['cL1']*p['cc7']*p['w0']*y[2]*(1 - y[5]) - p['cc5']*p['y_0']*(1-y[6])*y[5]

    #wr6 = p['cc32']*p['o0']*(1 - y[6])*y[3] + p['cc5']*(1 - y[6])*y[5] - (p['cc01']/p['y_0'])/tau(y[6],y[7],y[9],p=p)
    wr6 = p['cc32'] * p['o0'] * (1 - y[6]) * y[3] + p['cc5'] * (1 - y[6]) * z[1] - (p['cc01'] / p['y_0']) / tau(y[6],y[7],y[9],p=p) # Detailed PSI

    wr7 = 0.5*p['cc31']*p['o0']*(1-y[7])*y[3] + 0.5*(p['cc01']/p['z0'])/tau(y[6],y[7],y[9],p=p) - 0.5*y[8]*y[7]*(2*p['cc33']*p['n0']*y[11] + p['cL2']*Light2(y[9], p)*p['cc10']*p['q0']*(1 - y[4]) + p['cc8']*p['w0']*(1-y[2])) - p['cc36']*p['ccsu']*y[7]

    wr8 = p['cc11'] * psi2(y[8], y[9], rho, p) - p['cc11s'] * psis(y[8], p) + p['cc11h'] * p['s0'] * (1 - y[1]) * psih(y[8], y[9], y[12], y[1], p, rho) \
    - p['z0'] * y[8] * y[7] * (p['cc33'] * p['n0'] * y[11] + p['cc8'] * p['w0'] * (1 - y[2]) + L2 * Light2(y[9], p) * p['cc10'] * p['q0'] * (1 - y[4])) + p['cc1r'] * p['r0'] * y[8] * (1 - y[12]) \
    - p['cc19'] * p['n0'] * y[8] * (1 - y[11] - y[10]) + p['cc18'] * p['n0'] * y[11] - p['cc6'] * p['o0'] * p['w0'] * y[3] * (1 - y[2]) - 2 * p['cc1h'] * p['n0'] * p['s0'] * (1 + psig(y[3], p)) * psicc(y[11], y[1], y[8], y[3], y[12], p) \
    - 0 * 0.5 * p['cc17'] * p['o0'] * p['n0'] * y[3] * y[11] - p['cc31'] * p['o0'] * p['z0'] * y[3] * (1 - y[7]) - p['cc32'] * p['o0'] * p['y_0'] * y[3] * (1 - y[6])

    wr9 = p['cc30'] * p['q0'] * y[4] + p['cc15'] * p['cc01'] / tau(y[6], y[7], y[9], p=p) + p['cc31'] * p['o0'] * p['z0'] * y[3] * (1 - y[7]) - p['cc11'] * psi2(y[8], y[9], p=p, rho=rho) - p['cc11h'] * p['s0'] * (1 - y[1]) * psih(y[8], y[9], y[12], y[1], rho=rho, p=p)

    wr10 = p['cc1h']*p['s0']*(1+p['ccb']*psig(y[3],p=p))*psicc(y[11],y[1],y[8],y[3],y[12],p=p) +p['cc33']*p['z0']*y[8]*y[7]*y[11] + 0.5*p['cc17']*p['o0']*y[3]*y[11] - y[10]*(0.5*p['cc34']*p['w0']*(1-y[2]) + p['cc35'])

    wr11= p['cc19']*y[8]*(1 - y[11] - y[10]) - y[11]*(p['cc18'] + p['cc33']*p['z0']*y[8]*y[7] + 0.5*p['cc17']*p['o0']*y[3]) - p['cc1h']*p['s0']*(1+p['ccb']*psig(y[3],p=p))*psicc(y[11],y[1],y[8],y[3],y[12],p=p)

    wr12= p['cc1r']*y[8]*p['w0']*(1-y[2])*(1 - y[12]) - p['cc2r']*y[12]/(p['cc5r'] + p['cc6r']*p['r0']*(1 - y[12])) - p['cc3r']*y[12]/(p['cc7r']+p['cc8r']*p['r0']*(1 - y[12])) - p['cc4r']*p['o0']*y[3]*y[12]

    dydt1 = [wr1, wr2, wr3, wr4, wr5, wr6, wr7, wr8, wr9, wr10, wr11, wr12]

    # 1      2      3      4       5       6       7       8       9

    # x0    x1    1-a0    1-b0    1-a1    1-b1    1-fx    1-f      mv

    vr1 = p['cL1'] * p['cxx'] * (1 - z[1] - z[2]) - p['txx'] * z[2] - z[1] * (p['tx1'] * z[3] + p['tx2'] * z[5] + p['tx3'] * z[4] + p['tx4'] * z[6] + p['tx5'] * z[7] + p['tx6'] * z[8] + p['cc5'] * (1 - y[6]))

    vr2 = p['cL1'] * p['cxx'] * (1 - z[1] - z[2]) - z[2] * (p['txx'] + p['ca0'] * (1 - z[3]) + p['cb0'] * (1 - z[4])) + z[1] * (p['ta0'] * z[3] + p['tb0'] * z[4])

    vr3 = - z[3] * (p['ta0'] * z[1] + p['tx1'] * z[1] + p['ca1'] * (1 - z[5])) + (1 - z[3]) * (p['ca0'] * z[2] + p['ta1'] * z[5])

    vr4 = - z[4] * (p['tb0'] * z[1] + p['tx3'] * z[1] + p['cb1'] * (1 - z[6])) + (1 - z[4]) * (p['cb0'] * z[2] + p['tb1'] * (1 - z[6]))

    vr5 = - z[5] * (p['ta1'] * (1 - z[3]) + p['cfa'] * (1 - z[7]) + p['tx2'] * (z[1]) + p['caa'] * z[9]) + (1 - z[5]) * (p['ca1'] * z[3] + p['tfa'] * z[7])

    vr6 = - z[6] * (p['tb1'] * (1 - z[3]) + p['cfb'] * (1 - z[7]) + p['tx4'] * z[1] + p['cbb'] * z[9]) + (1 - z[6]) * (p['cb1'] * z[4] + p['tfb'] * z[7])

    vr7 = - z[7] * (p['tfa'] * (1 - z[5]) + p['tfb'] * (1 - z[6]) + p['cff'] * (1 - z[8]) + p['tx5'] * z[1]) + (1 - z[7]) * (p['cfa'] * z[5] + p['cfb'] * z[6] + p['tff'] * z[8])

    vr8 = - z[8] * (p['tff'] * (1 - z[7]) + p['tx6'] * z[1] + p['cmv'] * z[9] + p['cc7'] * y[2]) + (1 - z[8]) * (p['cff'] * z[7] + p['tmv'] * (1 - z[9]))

    vr9 = - z[9] * (p['cmv'] * z[8] + p['caa'] * z[5] + p['cbb'] * z[6]) + (1 - z[9]) * (p['tmv'] * (1 - z[8]) + p['coo'] * y[3])

    dydt2 = [vr1, vr2, vr3, vr4, vr5, vr6, vr7, vr8, vr9]

    return equinit(dydt1+dydt2, p)


def main():
    print('KEK')


if __name__ == '__main__':
    print('Hi')
