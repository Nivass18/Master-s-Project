#
# Truss Optimization with Genetic Algorihm
#

import numpy as np
import math
from numpy import random
from scipy.optimize import minimize
import sys
import time

try:
    population_size=int(input('Enter the population size ')) # Population size from the user
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)
try:
    mutation_probability = float(input('Enter the probability of mutation ')) # Mutation probability for Genetic Algorithm
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)

Target_value=25

x=np.zeros(population_size)                     # Initial Array for X values 
y=np.zeros(population_size)                     # Initial Array for Y values

init_population=np.zeros((population_size,2))   # Initial array for the initial population 
weights=np.zeros((population_size))             # Initial Weight arrays
str_lt=np.zeros((population_size))              # Initial limits array 
children=np.zeros((population_size,2))          # Initial array for the children

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

IP=1                                                                    # Value for Genetic Algorithm penalty function 
P2=1000                                                                 # Value for Genetic Algorithm penalty function

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def initial_population(population_size):
    '''
    ======================================================================
    Initial Population Function
    ----------------------------------------------------------------------
    population_size: Total number of population
    ----------------------------------------------------------------------
    return: 
    Init_population: Random array of parameters 
    ======================================================================
    '''
    #
    # Random values of X and Y values using random number generator
    #
    for i in range(population_size):
            x[i]=round(np.random.uniform(low=0, high=5),1)
            y[i]=round(np.random.uniform(low=0, high=5),1)
            init_population[i]=[x[i],y[i]]
    #
    # End
    #
    return init_population
#-----------------------------------------------------------------------------------------------------
#Calculating Fitness values---------------------------------------------------------------------------
def fitness_function(init_population):
    '''
    ======================================================================
    Fitness calculation for each gene in Population
    ----------------------------------------------------------------------
    init_population: Inital Population
    ----------------------------------------------------------------------
    return: 
    pf: Penalty function
    str_lt: Constraint value of each population
    ======================================================================
    '''
    #
    # Weight calculation of each gene in the population and the penalty function
    # 
    for i in range(population_size):
        weights[i]=((init_population[i][0]*init_population[i][1]))
        str_lt[i]=((2*init_population[i][0]+2*init_population[i][1])-20)
        pf_const=P2*(str_lt)
        pf=weights+pf_const
    #
    # End
    #

    return pf,str_lt
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def scaled_fitness(pf_value_array):
    '''
    ======================================================================
    Scaled fitness calculation function
    ----------------------------------------------------------------------
    pf_value_array: Penalty function value
    ----------------------------------------------------------------------
    return: 
    pf_array_sorted: Penalty function sorted array
    pf_array_sorted_rank: Rank of pealty function with sorted
    ======================================================================
    '''
    pf_array_sorted=np.sort(pf_value_array)[::-1]                       # Sorting the penalty array
    pf_array_sorted_elements=np.arange(1,(np.size(pf_array_sorted)+1))
    rank_sum=np.sum(pf_array_sorted_elements)                           # Calculating the rank of penalty function array
    pf_array_sorted_rank=pf_array_sorted_elements/rank_sum
    return pf_array_sorted,pf_array_sorted_rank
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def roulette_wheel_pop(input_array, pf_array_sorted_rank, number):
    '''
    ======================================================================
    Parent Seletion ( Roulette wheel )
    ----------------------------------------------------------------------
    input_array: Input array of population
    pf_array_sorted_rank: Rank of pealty function with sorted
    number: Number of parent to be seleted for Genetic Algorithm
    ----------------------------------------------------------------------
    return: 
    chosen_array: Array of Selected parents
    ======================================================================
    '''
    chosen = []                          # Initial Array for selecting parent
    max_value=max(pf_array_sorted_rank)
    #
    # Randomly select parent baed on the random value
    #
    for n in range(number):
        r = np.random.uniform(0,max_value)
        #
        # Selecing parent based on the random value in the population
        #
        for (i, individual) in enumerate(input_array):
            if r <= pf_array_sorted_rank[i]:
                chosen.append((individual))
                break
        #
        # End
        #
        chosen_array=np.array(chosen)
    #
    # End
    #
    return chosen_array
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def single_pt_crossover(A,B,x):
    '''
    ======================================================================
    Crossover function
    ----------------------------------------------------------------------
    A: Parent 1
    B: Parent 2
    x: Crossover point
    ----------------------------------------------------------------------
    return: 
    A_new: Children 1
    B_new: Children 2
    ======================================================================
    '''
    A_new=np.zeros(np.shape(A))     # Initial array for children_1
    B_new=np.zeros(np.shape(B))     # Initial array for children_2
    # Crossover operation based on single point 
    A_new=np.append(A[:x],B[x:])
    B_new=np.append(B[:x],A[x:])
    return A_new,B_new
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def multi_pt_crossover(A,B,X):
    '''
    ======================================================================
    Crossover function
    ----------------------------------------------------------------------
    A: Parent 1
    B: Parent 2
    X: Array of Crossover point
    ----------------------------------------------------------------------
    return: 
    A: Children 1
    B: Children 2
    ======================================================================
    '''
    #
    # Crossover operation for selected two parents and multi points
    #
    for i in X:
        A,B=single_pt_crossover(A,B,i)  # function call for single point crossover
    #
    # End
    #
    return A,B
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def randomly_mutate_population(population, mutation_probability):
    '''
    ======================================================================
    Mutation Function
    ----------------------------------------------------------------------
    population: Total population
    mutation_probability: probability of mutation
    ----------------------------------------------------------------------
    return: 
    population: Mutated population
    ======================================================================
    '''
    random_mutation_array = np.random.uniform(size=(population.shape))      # Create Random number
    random_mutation_boolean = random_mutation_array <= mutation_probability # Convert the random number to boolean operation
    
    #
    # Randomly mutate the respective selected parameter and update the value within the limits
    #
    for i in range(len(population)):
        if(random_mutation_boolean[i][0]):
            population[i][0]=round(np.random.uniform(low=0, high=5),1)
        if(random_mutation_boolean[i][1]):
            population[i][1]=round(np.random.uniform(low=0, high=5),1)
    #
    # End
    #
    return population
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def GA(population,pf):
    '''
    ======================================================================
    Genetic Algorithm
    ----------------------------------------------------------------------
    population: Total population
    pf: Penalty function of each gene in the population
    ----------------------------------------------------------------------
    return: 
    mutated_population: Population with GA operations applied
    ======================================================================
    '''
    #
    # Parent selection based on Roulette Wheel Selection
    #
    for k in range(int(population_size/2)):
        parent=roulette_wheel_pop(population, pf_array_sorted_rank, number=2)

        parent_1=parent[0]
        parent_2=parent[1]


        X=[]    # Initial array for creating crossover points
        #
        # Randomly create crossover points for the crossover operation
        #
        for i in range(2):
            num=np.random.randint(0,len(parent_1))
            if num not in X: 
                X.append(num)
                X=sorted(X)
        X=np.array(X)
        #
        # End
        #

        A_new,B_new=multi_pt_crossover(parent_1,parent_2,X) # Multipoint crossover function call
        children[2*k]=A_new
        children[(2*k)+1]=B_new

    
        mutated_population = randomly_mutate_population(children, mutation_probability) # Mutation function call 

        #
        # Eltisim and replace best fittest gene in the population
        #
        for i,x in enumerate(pf):
            for j in range(int(population_size/2)):
                if x==pf_array_sorted[j]:
                    mutated_population[0]=population[i]
        #
        # End
        #
    #
    # End
    #
    return mutated_population
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
Points=[]       # Initial array for getting the Penalty function value of each gene
def print_population(population, generation_number):
    '''
    ======================================================================
    Printing Results
    ----------------------------------------------------------------------
    population: Total population
    generation_number: Number of generations
    ----------------------------------------------------------------------
    ======================================================================
    '''
    print('\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
    print("Generation #",generation_number, ' Fittest Chromosome Fitness:',pf_array_sorted[0])
    print('Target Value', Target_value)
    print('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
    i=0
    for x in population:
        print("Chromosome #",i, ' :',x, '| Fitness: ',pf_array_sorted[i])
        print('---------------------------------------------------------')
        i=i+1
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

population=initial_population(population_size)  # Population generation

print('Initial Population \n',population)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
start = time.time()
pf,str_lt=fitness_function(population)  # Getting the fittness values and Penalized fitness based on constraints

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
pf_array_sorted,pf_array_sorted_rank=scaled_fitness(pf) # Sorting the rank and sorting the penalized fitness
population_1=np.zeros(np.shape(population))             # Initial array for next population generation
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# Rearranging the population based on the fitness
#
for i,x in enumerate(pf_array_sorted):
    for j in range(len(pf)):
        if pf[j]==x:
            population_1[i]=population[j]
#
# End
#
print_population(population_1,0) # Printing the population

generation_number=1

Points.append(pf_array_sorted[0])   # Append the penalty function value
#
# Loop for genetic algorithm parent selection, crossover and mutation
#
while pf_array_sorted[0]-Target_value:
    population=GA(population_1,pf)                          # GA algorithm function call
    pf,str_lt=fitness_function(population)                  # Getting the fittness values and Penalized fitness based on constraints
    pf_array_sorted,pf_array_sorted_rank=scaled_fitness(pf) # Sorting the rank and sorting the penalized fitness
    #
    # Rearranging the population based on the fitness
    #
    for i,x in enumerate(pf_array_sorted):
        for j in range(len(pf)):
            if pf[j]==x:
                population_1[i]=population[j]
    #
    # End
    #
    print_population(population_1,generation_number)    # Printing the each population

    Points.append(pf_array_sorted[0])
    generation_number+=1


print('The best solution is \n', population_1[0])
end = time.time()
temp = end-start
hours = temp//3600
temp = temp - 3600*hours
minutes = temp//60
seconds = temp - 60*minutes
print('Total time elapsed: %d:%d:%d' %(hours,minutes,seconds))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def objective(x, sign=-1.0):
    '''
    ======================================================================
    Objective function
    ----------------------------------------------------------------------
    x: Values of x
    ----------------------------------------------------------------------
    return: 
    X1*X2
    ======================================================================
    '''
    x1 = x[0]
    x2 = x[1]
    return sign*(x1*x2)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def constraint1(x):
    '''
    ======================================================================
    Constraints function
    ----------------------------------------------------------------------
    x: Values of x
    ----------------------------------------------------------------------
    return: 
    constrain values 
    ======================================================================
    '''
    return 20.0 - 2*x[0] -2*x[1]

x0 = [0,0]

b1 = (1,5)
b2 = (1,5)
bnds= (b1,b2)
con1 = {'type': 'ineq', 'fun': constraint1}
cons = [con1]
sol = minimize(objective, x0, method='SLSQP', bounds = bnds, constraints = cons)    # Function call to optimize the problem with constraints

print(sol)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
with open('Genetic_algorithm_with_Population_size_100_0.7.txt','w') as filename:     # Writing the results in a file for plotting
    for item in Points:
        filename.write('%s\n' %item)
