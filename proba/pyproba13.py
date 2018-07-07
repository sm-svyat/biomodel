#!/usr/bin/python3
import time
import json
from equation_parts import psih, psig, psicc, psi2, psis, Light2, tau, alpha, beta
from math import log10, exp
from scipy.integrate import ode
from engines import ResultsHandler, Step


def read_constants(const_path):
    '''
    Функция читает константы из ini-файла матлабовской модели
    :param CONST_PATH:
    :return:
    '''
    filename_extension = const_path.split('.')[1]
    constants = dict()
    with open(const_path, 'r') as const_file:
        if filename_extension == 'txt':
            for row in const_file:
                try:
                    const_name, const_values = row.split('#')
                except ValueError:
                    continue
                const_value = const_values.split('%')
                constants[const_name.split(' ')[0]]=float(const_value[0])
        elif filename_extension == 'json':
            json_constants = json.loads(const_file.read())
            constants = {key: json_constants[key]['value'] for key in json_constants}

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


def light_mode_switch(params, time, dataContainer):
    '''
    Функция переключения освещения. Константа мощности освещения меняется в зависимости от значений функции хевисайда
    в текущий момент времени. Функция хевисайда отражает режим освещения.
    :param params:
    :param time:
    :param dataContainer:
    :return new params:
    '''
    lightONParametrs = dataContainer.lightinModeDecloration['LightONParametrs']
    for key in lightONParametrs:
        params[key] = lightONParametrs[key]*dataContainer.heaviside()(time)
        #print(params[key])


    return params


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


#def equinit(y, p, equ_switch=[1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], equZ_switch=[0, 0, 0, 0, 0, 0, 0, 0, 0]):
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


def runge_kutta_method(t, y, fx, hs, p):
    '''
    Classical Runge–Kutta method (RK4)
    :param t:
    :param y:
    :param fx:
    :param hs:
    :return [y_val + (k1_val + 2*(k2_val + k3_val) + k4_val)/6 for y_val, k1_val, k2_val, k3_val, k4_val in zip(y,k1,k2,k3,k4)]:
    '''
    fx_res = fx(t, y, p)
    k1 = [res * hs for res in fx_res]
    yk = [y_val + k1_val*0.5 for y_val, k1_val in zip(y, k1)]
    fx_res2 = fx(t, yk, p)
    k2 = [res * hs for res in fx_res2]
    yk = [yk_val + k2_val*0.5 for yk_val, k2_val in zip(yk, k2)]
    fx_res3 = fx(t, yk, p)
    k3 = [res * hs for res in fx_res3]
    yk = [yk_val + k3_val for yk_val, k3_val in zip(yk, k3)]
    fx_res4 = fx(t, yk, p)
    k4 = [res * hs for res in fx_res4]
    return [y_val + (k1_val + 2*(k2_val + k3_val) + k4_val)/6 for y_val, k1_val, k2_val, k3_val, k4_val in zip(y,k1,k2,k3,k4)]


def rk4_method(y, t0, dataContainer, tf, p):
    '''
    Функция запуска метода Рунге-Кутты четвертого прядка
    Реализован данный метод в функцие runge_kutta_method
    :param y:
    :param t0:
    :param x:
    :param time:
    :param tf:
    :return x, time:
    '''
    t = t0
    step = Step(tStart=t0, tEnd=tf)

    if dataContainer.lightinModeDecloration:
        while t < tf:  # Need coroutine maybe FIXME
            p = light_mode_switch(p, t, dataContainer)

            try:
                y = runge_kutta_method(t, y, f, step.step, p)
            except OverflowError:
                print("Rise OverflowError: math range error, when t = %s" % t)
                break

            boundary_conditions(y)

            # Сохранение данных полученных данных в контейнер
            try:
                dataContainer.timeNpList[dataContainer.stepGrid.currentStep] = t
                for i in range(len(dataContainer.variablesList)):
                    dataContainer.variablesNpList[i, dataContainer.stepGrid.currentStep] = y[i]
            except IndexError:
                print("IndexError")
            dataContainer.stepGrid.currentStep += 1
            t += dataContainer.stepGrid.step

    else:
        while t < tf: #Need coroutine maybe FIXME
            try:
                y = runge_kutta_method(t, y, f, step.step, p)
            except OverflowError:
                print("Rise OverflowError: math range error, when t = %s" %t)
                break
            except ValueError:
                print('Уменьние шага в N раз.\nN = {}'.format(step.reductionFactor))
                 #step.reduce()
                continue
            boundary_conditions(y)

            try:
                dataContainer.timeNpList[dataContainer.stepGrid.currentStep] = t
                for i in range(len(dataContainer.variablesList)):
                    dataContainer.variablesNpList[i, dataContainer.stepGrid.currentStep] = y[i]
            except IndexError:
                print("IndexError")
            dataContainer.stepGrid.currentStep += 1
            t += dataContainer.stepGrid.step
    return dataContainer


def ode15s_solver(y, t0, dataContainer, tf, p, ode_es):
    '''
    Функция запуска ode решателя ode15s из библиотеки scipy
    :param y:
    :param t0:
    :param x:
    :param time:
    :param tf:
    :return x, time:
    '''
    t = t0
    step = Step(tStart=t0, tEnd=tf)
    #ode15s = ode(f)
    ode15s = ode(ode_es)
    #ode15s.set_integrator('vode', method='bdf', order=5, nsteps=3000, atol=1e-4, rtol=1e-4)
    #ode15s.set_integrator('vode', method='adams')
    ode15s.set_integrator('lsoda')
    ode15s.set_f_params(p)
    ode15s.set_initial_value(y, t0)

    #Проверка режима освещения
    if dataContainer.lightinModeDecloration:
        while ode15s.successful() and ode15s.t < tf:
            p = light_mode_switch(p, ode15s.t, dataContainer)
            #ode15s.set_f_params(p)
            try:
               y = ode15s.integrate(ode15s.t+step.step)
            except OverflowError:
                print("Rise OverflowError: math range error, when t = %s" %t) #
                break
            except ValueError:
                print('Уменьние шага в N раз.\nN = {}'.format(step.reductionFactor))
                step.reduce()
                continue
            boundary_conditions(y)
            try:
                dataContainer.timeNpList[dataContainer.stepGrid.currentStep] = ode15s.t
                for i in range(len(dataContainer.variablesList)):
                    dataContainer.variablesNpList[i, dataContainer.stepGrid.currentStep] = y[i]
            except IndexError:
                print("IndexError")
            dataContainer.stepGrid.currentStep += 1

    else:
        while ode15s.successful() and ode15s.t < tf:
            try:
                y = ode15s.integrate(ode15s.t + step.step)
            except OverflowError:
                print("Rise OverflowError: math range error, when t = %s" % t)  #
                break
            except ValueError:
                print('Уменьние шага в N раз.\nN = {}'.format(step.reductionFactor))
                #step.reduce()
                continue
            boundary_conditions(y)
            # dataContainer.add_variables(ode15s.t, y)

            try:
                dataContainer.timeNpList[dataContainer.stepGrid.currentStep] = ode15s.t
                for i in range(len(dataContainer.variablesList)):
                    dataContainer.variablesNpList[i, dataContainer.stepGrid.currentStep] = y[i]
            except IndexError:
                print("IndexError")
            dataContainer.stepGrid.currentStep += 1
    return dataContainer


def adaptiv_step_solver(y, t, dataContainer, tf):
    '''
    Функция запуска ode решателя ode45 из библиотеки scipy
    Использует переменный шаг
    Нуждается в доработке
    :param y:
    :param t:
    :param x:
    :param time:
    :param tf:
    :return x, time:
    '''
    backend = 'dopri5'
    solver = ode(f).set_integrator(backend, nsteps=1000, atol=1e-4, rtol=1e-4)

    def solout(t, y):
        boundary_conditions(y)
        dataContainer.add_variables(t, y)

    solver.set_solout(solout) #при каждом вычислении выполняет функцию solout
    solver.set_initial_value(y, t)
    solver.integrate(tf)
    return dataContainer


def adaptiv_step_solver2(y, t, data_container, tf):
    # backend = 'dopri5'
    backend = 'dop853'
    solver = ode(f).set_integrator(backend, nsteps=3000, atol=1e-8, rtol=1e-8)

    def solout(t, y):
        data_container.add_variables(t, y)
        print('hi')

    solver.set_solout(solout)  # при каждом вычислении выполняет функцию solout
    solver.set_initial_value(y, t)

    solver.integrate(tf)
    return y


if __name__ == '__main__':
    startTime = time.time() # нужно для расчета времени вычислений
    #Определение констант и начальных значений переменных
    CONST_PATH = 'Proba13.ini'
    p = read_constants(CONST_PATH)

    #Переопределение некоторых констант
    p = equcoeff(p)
    t0 = 0 #p['t0'] # Начальное время
    tf = 8000 #p['tend'] # Конечное время

    #Создание списка наименований переменных
    y0names = ['y0(' + str(i) + ')' for i in range(1, 13)]
    y0 = [p[name] for name in y0names]
    z0names = ['zz0(' + str(i) + ')' for i in range(1, 10)]
    z0 = [p[name] for name in z0names]
    y = y0+z0

    #Создание класса для записи промежуточных результатов. Можно предварительно задать режим освещения в lightinModeDecloration, если требуется.
    dataContainer = ResultsHandler(t0, tf, y)
    #dataContainer = ResultsHandler(t0, y, lightinModeDecloration = {'TimeLightON': [0, 300, 600, 900], 'TimeLightOFF': [2, 302, 602, 902], 'LightONParametrs': {'cL1': 1}, 'LightOFFParametrs': {'cL1': 0}})
    #dataContainer = ResultsHandler(t0, tf, y, lightinModeDecloration = {'TimeLightON': [0], 'TimeLightOFF': [12000], 'LightONParametrs': {'cL1': 0.006}, 'LightOFFParametrs': {'cL1': 0}})

    #Создание выходного файла
    dataContainer.create_new_output()

    #Добавление параметров в класс с промежуточными результатами
    dataContainer.parameters = p

    #Словарь доступных ode решателей. Работают ode15s и rkm. ode15s быстрее.
    solvers = {'rkm': rk4_method, 'adaptiv_step': adaptiv_step_solver, 'ode15s': ode15s_solver}

    #Запуск определенного ode решателя
    readyCalculation = solvers['ode15s'](y, t0, dataContainer, tf, p)

    #Сохранение реузультатов в файл
    readyCalculation.save_results()

    endTime = time.time() # нужно для расчета времени вычислений
    print('Время выполнения кода {} минут(ы)'.format(round((endTime - startTime)/60, 2)))

    #Вывод интересуещих кривых
    #readyCalculation.graph(["x0"], exp_data_flag=['llvw'], parametrs_flag=True)
    #readyCalculation.graph(["x0", "z"], exp_data_flag=['llvw'])
    readyCalculation.graph(["x0"], parametrs_flag=True)
