import sys
import pandas as pd
import matplotlib.pyplot as plt

def graph(data_frame, x, plot_variable=['x0'], exp_data_flag=False, parametrs_flag=False, plot_graph=False, save_graph_as=False):
    '''                                                                          
    Функция строит графики для переменных указанных в plot_variable              
    :param plot_variable: список переменных                                      
    :param exp_data_flag: путь к файлу с экспериментальными данными              
    :param parametrs_flag: если True выводит на графике значения парметров модели
    :return:                                                                     
    '''                                                                          
    plt.figure(figsize=(10, 6))                                                  
    for variable in plot_variable:                                               
        curve = plt.plot(data_frame[x], data_frame[variable],'-', alpha=0.8, label=variable)
        plt.setp(curve, linewidth=3)                                             
        plt.tick_params(axis='both', labelsize=20)                               
                                                                                 
    plt.grid(True)                                                               
    plt.legend(loc='center', bbox_to_anchor=(1.05, 0.5), shadow=True)            
    if x == 't':
        plt.xlabel('$time, s$', size=15)                                             
    plt.title('Graph of the relative carriers concentrations.\nModel of detailed PSI', size=16)
                                                                                 
    if plot_graph:                                                               
        plt.show()                                                               
                                                                                 
    if save_graph_as:                                                            
        plt.savefig(save_graph_as)                                      


def help():
    print('python plot.py <file_name> <graph_name> <plot_variables y1,y2,y3>')


if __name__ == '__main__':
    try:
        file_name, graph_name, plot_variables = sys.argv[1:4]
    except ValueError:
        print('Передано некорректное количество аргументов')
        help()
        raise SystemExit(1)

    #file_name = 'results/my_result.txt'
    try:
        df = pd.read_csv(file_name, sep='\t', skiprows=range(1,2))
        #data = read_data(file_name)
    except FileNotFoundError:
        print('Не удалось открыть файл {}'.format(file_name))
        help()
        raise SystemExit(1)
    
    plot_variables = plot_variables.split(',')

    try:
        graph(df, plot_variables[0], plot_variables[1:], save_graph_as=graph_name)
        #data_storage.graph([variable], save_graph_as=graph_name)
    except KeyError:
        print('Указанно некорректное название переменной {}'.format(variable))

