import copy
from math import log10
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


class ResultsHandler:
    def __init__(self, initialTime, tEnd, initialValues, lightinModeDecloration=None):
        self.initialValues = initialValues
        self.initialTime = initialTime
        self.variablesList = [[variable] for variable in initialValues]
        self.timeList = [initialTime]
        self.headers = ["s","w","o","q","x","y","z","u","v","n","m","r", "x0", "x1", "1-a0", "1-b0", "1-a1", "1-b1", "1-fx", "1-f", "mv"]
        self.variables_dict = {key: value for key, value in zip(self.headers, range(len(self.headers)))}
        self.lightinModeDecloration = lightinModeDecloration
        self.parameters = None
        self.data_initialization(initialTime, tEnd)

    def data_initialization(self, initialTime, tEnd):
        self.stepGrid = Step(tStart=initialTime, tEnd=tEnd)
        self.nsteps = tEnd / self.stepGrid.step
        self.timeNpList = np.zeros(int(round(self.nsteps, 0)))
        self.variablesNpList = np.zeros((len(self.initialValues), int(round(self.nsteps, 0))))
        self.stepGrid.currentStep = 0

    def __repr__(self):
        return str(self.variablesNpList)

    def add_variables(self, time, variables):
        '''
        Функция добавления новых значений переменных
        :param time: текущее время
        :param variables: текущие значения переменных
        :return:
        '''
        try:
            self.timeNpList[self.stepGrid.currentStep] = time
            for i in range(len(self.variablesList)):
                self.variablesNpList[i, self.stepGrid.currentStep] = variables[i]
        except IndexError:
            print("IndexError")
        self.stepGrid.currentStep += 1

    def create_new_output(self, filename = 'my_out.txt'):
        '''
        Функция создает файл для записи посчитанных значений. Создает в файле строку переменных и строку начальных значений.
        :param filename: Имя файла. default: 'my_out.txt'
        :return:
        '''
        with open(filename, 'w', encoding='utf-8') as output:  # creating new out file with headers all ready
            output.write('t\t' + '\t'.join(self.headers) + '\n')
            output.write('\t'.join(
                ['time', '[ATP]', '[Fd]', '[O2]', '[PS2]', '[PS1]', '[Pc]', '[QH2]', '[pHout]', '[pHin]', '[NADP]',
                 '[NADPH+]', '[Tr]', "[P700+]", "[P700*]", "[1-a0]", "[1-b0]", "[1-a1]", "[1-b1]", "[1-fx]", "[1-f]", "[mv]"]) + '\n')
        with open(filename, 'a', encoding='utf-8') as output:
            ready_results = copy.deepcopy(self.initialValues)
            for i in [7, 8]:
                ready_results[i] = 6 - log10(ready_results[i])
            self.initialStringValues = [str(i) for i in ready_results]
            output.write(str(self.initialTime / 100) + "\t" + '\t'.join(self.initialStringValues) + '\n')

    def save_results(self, filename='my_out.txt'):
        '''
        Функция записывает текущий набор данных в файл
        :param filename: Имя файла. default: 'my_out.txt'
        :return:
        '''
        self.create_new_output(filename)
        ready_results = copy.deepcopy(self.variablesNpList)
        scaledTimeList = [time / 100 for time in self.timeNpList]
        for i in [7, 8]:
            ready_results[i] = [6 - log10(variable) for variable in self.variablesNpList[i]]
        with open(filename, 'a', encoding='utf-8') as output:
            for point in [i for i in range(len(self.timeNpList)) if i%(self.nsteps/1000) == 0]:
                interim_results = [str(i[point]) for i in ready_results]
                output.write(str(scaledTimeList[point]) + "\t" + '\t'.join(interim_results) + '\n')


    def graph(self, plot_variable=['x'], exp_data_flag=False, parametrs_flag=False, plot_graph=False, save_graph_as=False):
        '''
        Функция строит графики для переменных указанных в plot_variable
        :param plot_variable: список переменных
        :param exp_data_flag: путь к файлу с экспериментальными данными
        :param parametrs_flag: если True выводит на графике значения парметров модели
        :return:
        '''
        plt.figure(figsize=(10, 6))
        scaledTimeList = [time / 100 for time in self.timeNpList]
        for variable in plot_variable:
            curve = plt.plot(scaledTimeList, self.variablesNpList[self.variables_dict[variable]], alpha=0.8, label=variable)
            plt.setp(curve, linewidth=3)    
            plt.tick_params(axis='both', labelsize=20)
        
        if exp_data_flag:
            for data in exp_data_flag:
                y_axis = list()
                x_axis = list()
                with open('../experimental data/' + data + '.txt', 'r') as exp_data:
                    for row in exp_data:
                        point = row.split(',')
                        x_axis.append(float(point[0]))
                        y_axis.append(float(point[1]))
                y_max = max(y_axis)
                y_axis = [y / y_max for y in y_axis]
                x_axis = [x for x in x_axis]
                #x_axis = [x/1000 for x in x_axis]
                plt.plot(x_axis, y_axis, alpha=0.8, label=data)
        plt.title('Graph of the relative carriers concentrations.\nModel of detailed PSI', size=16)

        if parametrs_flag:
            plt.text(self.stepGrid.tEnd*0.0071, 0.21, 'cc5 = {}\ncc7 = {}\n\nca0 = {}\ncb0 = {}\nca1 = {}\ncb1 = {}\ncfa = {}\ncfb = {}\ncff = {}\ncxx = {}\n\n'
                             'tx1 = {}\ntx2 = {}\ntx3 = {}\ntx4 = {}\ntx5 = {}\ntx6 = {}'.format(self.parameters['cc5'],
                                                                                                 self.parameters['cc7'],
                                                                                                 self.parameters['ca0'],
                                                                                                 self.parameters['cb0'],
                                                                                                 self.parameters['ca1'],
                                                                                                 self.parameters['cb1'],
                                                                                                 self.parameters['cfa'],
                                                                                                 self.parameters['cfb'],
                                                                                                 self.parameters['cff'],
                                                                                                 self.parameters['cxx'],
                                                                                                 self.parameters['tx1'],
                                                                                                 self.parameters['tx2'],
                                                                                                 self.parameters['tx3'],
                                                                                                 self.parameters['tx4'],
                                                                                                 self.parameters['tx5'],
                                                                                                 self.parameters['tx6'],))

        plt.grid(True)
        plt.legend(loc='center', bbox_to_anchor=(1.05, 0.5), shadow=True)
        plt.xlabel('$time, s$', size=15)
        
        if plot_graph: 
            plt.show()

        if save_graph_as:
            plt.savefig(save_graph_as)
        

    def concentrations(self):
        '''
        Функция строит pHout и pHin
        :return:
        '''
        plot_variable = ['u', 'v']
        plt.figure(figsize=(10, 5))
        scaledTimeList = [time / 100 for time in self.timeNpList]
        for variable in plot_variable:
            x = list(self.variablesNpList[self.variables_dict[variable]])
            for i in range(len(self.variablesNpList[self.variables_dict[variable]])):
                x[i] = 6 - log10(x[i])
            plt.plot(scaledTimeList, x, '-', alpha=0.8, label=variable)

        plt.legend(loc='center', bbox_to_anchor=(1.05, 0.5), shadow=True)
        plt.xlabel('$time, s$', size=14)
        plt.show()

    def heaviside(self):
        '''
        Генерирует функцию хевисайда согласно lightinModeDecloration
        :return:
        '''
        return lambda x: 0.5 * sum([np.sign(x - s) + 1 for s in self.lightinModeDecloration['TimeLightON']]
                                   + [-np.sign(x - s) - 1 for s in self.lightinModeDecloration['TimeLightOFF']])


    def graph_hs(self):
        '''
        Строит функцию хевисайда
        :return:
        '''
        plt.figure(figsize=(10, 6))
        plt.plot(self.timeList, [self.heaviside()(i) for i in self.timeList], alpha=0.8)
        plt.title('Graph of the lighting condition.', size=16)
        plt.xlabel('$time, 100s $', size=10)
        plt.show()



class Step:
    def __init__(self, tEnd, tStart=0, boostFactor=1, reductionFactor=1):
        self.tEnd = tEnd
        self.step = 0.01
        self.currentStep = 0
        self.boostFactor = boostFactor
        self.reductionFactor = reductionFactor

    def reduce(self):
        self.step = self.step/self.reductionFactor

    def increase(self):
        self.step = self.step*self.boostFactor
