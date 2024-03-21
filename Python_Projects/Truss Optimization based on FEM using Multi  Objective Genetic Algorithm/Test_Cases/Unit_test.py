import numpy as np
import sys
import math
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
    max_length_of_bar=int(input("Enter the length of the bar "))                                    # Maximum length of the bar

    number_of_span=int(input("Enter the number of span in the truss structure "))                   # Total number of spans in the truss 
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
nodes=np.zeros(2*(number_of_span+1))            # Initial array for nodes 
values=np.zeros((2,2*(number_of_span+1)))       # Initial array for node values 

x_coordinates=np.zeros((2*(number_of_span+1)))  # Initial array for X coordinate points 
y_coordinates=np.zeros((2*(number_of_span+1)))  # Initial array for Y coordinate points 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_all_elements_in_input_array(population):
    '''
    ======================================================================
    Get all the elements from the input array Function
    ----------------------------------------------------------------------
    population_size: Total number of population
    ----------------------------------------------------------------------
    return: 
    count: Total number of Elements in the truss
    ======================================================================
    '''
    count=0                         # Initial value for get the count of elements

    #
    # Counting the number of elements in the population
    #
    for element in population:
        count+=len(element)
    #
    # End
    #

    return count
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def dir_cosine(from_point, to_point):
    '''
    ======================================================================
    Directional cosine Function
    ----------------------------------------------------------------------
    from_point: Coordinate value of from node
    to_point: Coordinate value of to node
    ----------------------------------------------------------------------
    return: 
    x_cosine,y_cosine : Directional cosine of each element
    ======================================================================
    '''
    vector=math.sqrt((to_point[0]-from_point[0])**2+(to_point[1]-from_point[1])**2) # Denominator value 
    x_cosine=(to_point[0]-from_point[0])/vector # Directional cosine calculation
    y_cosine=(to_point[1]-from_point[1])/vector # Directional cosine calculation
    return(x_cosine,y_cosine)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def transformation_matrix(x_cosine,y_cosine):
    '''
    ======================================================================
    Transformation matrix function
    ----------------------------------------------------------------------
    x_cosine,y_cosine : Directional cosine of each element
    ----------------------------------------------------------------------
    return: 
    Transformation matrix of each element
    ======================================================================
    '''
    return np.array([[x_cosine,y_cosine,0,0],[0,0,x_cosine,y_cosine]])
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    A_new=np.zeros(np.shape(A)) # Initial array for children_1
    B_new=np.zeros(np.shape(B)) # Initial array for children_2
    #
    # Crossover operation based on single point 
    #
    for i in range(3):
        A_new[i]=np.append(A[i][:x],B[i][x:])
        B_new[i]=np.append(B[i][:x],A[i][x:])
    #
    # End
    #
    return A_new,B_new
#----------------------------------------------------------------------------------------------------------------------------------------------------------------


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
        A,B=single_pt_crossover(A,B,i) # function call for single point crossover
    #
    # End
    #
    return A,B
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
def force_vector_calc(number_of_node):
    '''
    ======================================================================
    Force Vector from the user
    ----------------------------------------------------------------------
    number_of_node: Total number of nodes in the Truss
    ----------------------------------------------------------------------
    return: 
    force: Vector of force in the truss
    ======================================================================
    '''
    force=[]    # Initial Force array 

    #
    # Getting the x and y values of force vector from the user and arranging in the force vector
    #
    for i in range(number_of_node):
        try:
            x_value=int(input("Enter the value of force vector in " + str(i) + " in X dir "))   # Value of Force vector in X direction
            y_value=int(input("Enter the value of force vector in " + str(i) + " in Y dir "))   # Value of Force vector in Y direction
        except ValueError:
            print("\n Only integer values are valid")
            sys.exit(0)
        
        force.append([x_value,y_value])                                                     # Append Force vector in Force array
    #
    # End
    #
    
    return force
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_restrained_dofs():
    '''
    ======================================================================
    Get restrained DOFs Function
    ----------------------------------------------------------------------
    return: 
    restrained_dofs_array: Restrained DOFs array
    ======================================================================
    '''
    try:
           size_of_restrained_dofs=int(input("Enter the size of the restrained dofs array "))  # Getting size of restrained dof from user
    except ValueError:
        print("\n Only integer values are valid")
        sys.exit(0)
    

    restrained_dofs_array=np.zeros((size_of_restrained_dofs))                           # Arranging the restrained dof with zeros

    #
    # Getting the restrained dof from the user and arranging
    #
    for i in range(size_of_restrained_dofs):
        try:
           restrained_dofs_array[i]=int(input("Enter the numbers of the restrained dofs in the truss "))
        except ValueError:
            print("\n Only integer values are valid")
            sys.exit(0)
    
    restrained_dofs_array=restrained_dofs_array.astype('int32')  
    #
    # End
    #

    return restrained_dofs_array
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
