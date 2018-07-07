import time
import importlib
import sys
sys.path.append('proba')
from proba.pyproba13 import read_constants, equcoeff, ode15s_solver
from proba.engines import ResultsHandler


def help():
    print('python calc.py <equation_path> <const_path> <output_file>\n')


if __name__ == '__main__':
    startTime = time.time() # нужно для расчета времени вычислений 
    
    try:
        equations_path, const_path, output_file = sys.argv[1:4]
    except ValueError:
        print('Передано некорректное количество аргументов\n')
        help()
        raise SystemExit(1)

    #Импорт системы дифференциальных уравнений
    package_name = equations_path.split('.')[0]
    package = '.'.join(package_name.split('/'))
    try:
        module = importlib.import_module(package)
        ode_equation_system = module.f
        #from . import f as ode_equation_system

    except ImportError:
        print('Ошибка загрузки модуля {}\n'.format(sys.argv[1]))
        help()
        raise SystemExit(1)

    #Определение констант и начальных значений переменных
    #CONST_PATH = 'proba/Proba13.ini'
    try:
        p = read_constants(const_path)                                               
    except FileNotFoundError:
        print('Ошибка выгрузки параметров из {}\n'.format(const_path))
        help()
        raise SystemExit(1)

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
                                                                                 
    #Создание выходного файла                                                    
    #dataContainer.create_new_output()                                            
                                                                                 
    #Добавление параметров в класс с промежуточными результатами                 
    dataContainer.parameters = p                                                 
                                                                                 
    #Запуск определенного ode решателя                                           
    readyCalculation = ode15s_solver(y, t0, dataContainer, tf, p, ode_es=ode_equation_system)            
                                                                                 
    #Сохранение реузультатов в файл                                              
    readyCalculation.save_results(output_file)
    
    endTime = time.time() # нужно для расчета времени вычислений                 
    print('Время выполнения кода {} минут(ы)'.format(round((endTime - startTime)/60, 2)))
                                                                                 
    #Вывод интересующих кривых                                                   
    #readyCalculation.graph(["x0"], exp_data_flag=['llvw'], parametrs_flag=True) 
    #readyCalculation.graph(["x0", "z"], exp_data_flag=['llvw'])                 
    #readyCalculation.graph(["x0"], parametrs_flag=True, save_graph_as='graphics/my_graph.png')
