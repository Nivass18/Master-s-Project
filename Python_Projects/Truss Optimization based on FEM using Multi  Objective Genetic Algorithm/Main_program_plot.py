'''
Fitness and Iteration number plot with Generation size 100 and varying mutation probability 0.6, 0.7, 0.8
'''

import matplotlib.pyplot as plt
import numpy as np
#-----------------------------------------------------------------------------------------------------
#Getting x values with same generationa and varying mutation probability------------------------------
x_generation_size_40_6 = np.loadtxt('Genetic_algorithm.txt', unpack=True)
#-----------------------------------------------------------------------------------------------------
#for counting the maximum number of lines in each file------------------------------------------------
file_6 = open("Genetic_algorithm.txt", "r")
#-----------------------------------------------------------------------------------------------------
#Line count for mutation probability 0.6--------------------------------------------------------------
line_count_6 = 0
for line in file_6:
    if line != "\n":
        line_count_6 += 1
file_6.close()
#-----------------------------------------------------------------------------------------------------
#creating number of iteration for each mutation probability-------------------------------------------
y_6=np.arange(line_count_6)
#-----------------------------------------------------------------------------------------------------
#Plotting Graph---------------------------------------------------------------------------------------

plt.title("Fittness and Iteration Plot")    
plt.plot(y_6,x_generation_size_40_6, c='b', label='rho_m=0.6')
plt.text(line_count_6,x_generation_size_40_6[-1],str(round(x_generation_size_40_6[-1],2)))
plt.legend()
   

plt.xlabel('Iteration')
plt.ylabel('Fitness Value')


plt.savefig('main_plot.jpg')