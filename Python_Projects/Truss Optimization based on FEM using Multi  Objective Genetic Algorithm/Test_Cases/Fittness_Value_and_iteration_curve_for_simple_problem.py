'''
Fitness and Iteration number plot with Generation size 100 and varying mutation probability 0.6, 0.7, 0.8
'''

import matplotlib.pyplot as plt
import numpy as np
#-----------------------------------------------------------------------------------------------------
#Getting x values with same generationa and varying mutation probability------------------------------
x_generation_size_40_6 = np.loadtxt('Genetic_algorithm_with_Population_size_40_0.6.txt', unpack=True)
x_generation_size_40_7 = np.loadtxt('Genetic_algorithm_with_Population_size_40_0.7.txt', unpack=True)
x_generation_size_40_8 = np.loadtxt('Genetic_algorithm_with_Population_size_40_0.8.txt', unpack=True)
#-----------------------------------------------------------------------------------------------------
#for counting the maximum number of lines in each file------------------------------------------------
file_6 = open("Genetic_algorithm_with_Population_size_40_0.6.txt", "r")
file_7 = open("Genetic_algorithm_with_Population_size_40_0.7.txt", "r")
file_8 = open("Genetic_algorithm_with_Population_size_40_0.8.txt", "r")
#-----------------------------------------------------------------------------------------------------
#Line count for mutation probability 0.6--------------------------------------------------------------
line_count_6 = 0
for line in file_6:
    if line != "\n":
        line_count_6 += 1
file_6.close()
#-----------------------------------------------------------------------------------------------------
#Line count for mutation probability 0.7--------------------------------------------------------------
line_count_7 = 0
for line in file_7:
    if line != "\n":
        line_count_7 += 1
file_7.close()
#-----------------------------------------------------------------------------------------------------
#Line count for mutation probability 0.8--------------------------------------------------------------
line_count_8 = 0
for line in file_8:
    if line != "\n":
        line_count_8 += 1
file_8.close()
#-----------------------------------------------------------------------------------------------------
#creating number of iteration for each mutation probability-------------------------------------------
y_6=np.arange(line_count_6)
y_7=np.arange(line_count_7)
y_8=np.arange(line_count_8)
#-----------------------------------------------------------------------------------------------------
#Plotting Graph---------------------------------------------------------------------------------------
fig,ax = plt.subplots(3)
ax[0].set_title('Fitness and Iteration plot with varying mutation probability')
ax[0].set_title("Fittness and Iteration Plot")    
ax[0].plot(y_6,x_generation_size_40_6, c='b', label='rho_m=0.6')
ax[0].text(line_count_6,x_generation_size_40_6[-1],str([line_count_6,x_generation_size_40_6[-1]]))
leg = ax[0].legend()
   
ax[1].plot(y_7,x_generation_size_40_7, c='b', label='rho_m=0.7')
ax[1].text(line_count_7,x_generation_size_40_7[-1],str([line_count_7,x_generation_size_40_7[-1]]))
leg = ax[1].legend()
    
ax[2].plot(y_8,x_generation_size_40_8, c='b', label='rho_m=0.8')
ax[2].text(line_count_8,x_generation_size_40_8[-1],str([line_count_8,x_generation_size_40_8[-1]]))
leg = ax[2].legend()

for i in ax.flat:
    i.set(xlabel='Iteration', ylabel='Fitness Value')


plt.savefig('Fitness_value_and_iteration_plot_for_Population_size_40_with_varying_mutation_probability.jpg')