#
# Truss Optimization with Genetic Algorihm
#

import numpy as np
import math 
import matplotlib.pyplot as plt
import sys
import time

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
    total_number_of_parameters=int(input("Enter the number of parameters in each bar "))            # Total number of parameters in each bar
    number_of_span=int(input("Enter the number of span in the truss structure "))                   # Total number of spans in the truss 

    max_length_of_bar=int(input("Enter the length of the bar "))                                    # Maximum length of the bar
    max_height_of_bar=int(input("Enter the height of the bar "))                                    # Maximum height of the bar

    A_min=float(input('Enter the min value of area '))                                                # Minimum value of the area
    A_max=float(input('Enter the max value of area '))                                                # Maximum value of the area

    population_size=int(input('Enter the population size '))                                        # Population size from the user

    value=int(input('Enter the number of constraints (Stress (1) / strain (2) or both (3)) '))      # User defined values for the Constraints
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)

try:
    ele_density=float(input("Enter the density of element "))                                 # Density of element
    ele_stiffness=float(input("Enter the stiffness of element "))                             # Stiffness of element
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)

if value==1:
    max_strain=0
    try:
        max_stress=int(input('Enter stress limit '))                                                # Allowable stress in the Truss (N/mm**2)
    except ValueError:
        print("\n Only integer values are valid")
        sys.exit(0)
    
elif value==2:
    max_stress=0
    try:
        max_strain=int(input('Enter strain limit '))                                                # Allowable strain in the Truss (mm)  
    except ValueError:
        print("\n Only integer values are valid")
        sys.exit(0)
else:
    try:
        max_stress=int(input('Enter stress limit '))                                                # Allowable stress in the Truss (N/mm**2)
        max_strain=int(input('Enter strain limit '))                                                # Allowable strain in the Truss (mm)
    except ValueError:
        print("\n Only integer values are valid")
        sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

IP=1                                                                    # Value for Genetic Algorithm penalty function 
P2=1000                                                                 # Value for Genetic Algorithm penalty function 
try:
    mutation_probability = float(input('Enter the probability of mutation ')) # Mutation probability for Genetic Algorithm
except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)
 

number_of_node=2*(number_of_span+1)                                     # Total number of nodes in the Truss

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

VP=np.zeros((population_size,number_of_span+1))     # Initial array for the Vertical Position of the each node
VL=np.zeros((population_size,number_of_span+1))     # Initial array for the Vertical Length of the each node
D=np.zeros((population_size,number_of_span+1))      # Initial array for the Diagonal bar of the Truss
AV1=np.zeros((population_size,number_of_span+1))    # Initial array for the Area of the Vertical Bar in the Truss
AB1=np.zeros((population_size,number_of_span+1))    # Initial array for the Area of the Bottom Bar in the Truss
AT1=np.zeros((population_size,number_of_span+1))    # Initial array for the Area of the Top Bar in the Truss
AD11=np.zeros((population_size,number_of_span+1))   # Initial array for the Area of the Bar going Downwards in the Truss
AD12=np.zeros((population_size,number_of_span+1))   # Initial array for the Area of the Bar going Upwards in the Truss

#----------------------------------------------------------------------------------------------------------------------------------------------------------------

Init_population=np.zeros((population_size,number_of_span+1,total_number_of_parameters)) # Initial array for the initial population 
children=np.zeros((population_size,number_of_span+1,total_number_of_parameters))        # Initial array for the children         

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def initial_population(max_length_of_bar,max_height_of_bar,number_of_span,A_min,A_max,population_size):
    '''
    ======================================================================
    Initial Population Function
    ----------------------------------------------------------------------
    max_length_of_bar: Maximum length of the bar in the truss
    max_height_of_bar: Maximum height of the bar in the truss
    number_of_span: Total number of span in the truss
    A_min: Minimum value of the cross-section area in each element of the truss
    A_max: Maximum value of the cross-section area in each element of the truss
    population_size: Total number of population
    ----------------------------------------------------------------------
    return: 
    Init_population: Random array of parameters 
    ======================================================================
    '''
    #
    # Create Random values for each parameter in the truss and return the initial population array
    #
    for i in range(population_size):
        for j in range(number_of_span+1):
            VP[i][j]=round(np.random.uniform(low=0, high=max_height_of_bar-max_length_of_bar))  # Random array for the Vertical Position of the each node
            VL[i][j]=round(np.random.uniform(low=max_length_of_bar, high=max_height_of_bar))    # Random array for the Vertical Length of the each node
            D[i][j]=round(np.random.randint(low=1.0, high=4.0))                                 # Random array for the Diagonal bar of the Truss
            AV1[i][j]=round(np.random.uniform(low=A_min, high=A_max))                           # Random array for the Area of the Vertical Bar in the Truss
            AB1[i][j]=round(np.random.uniform(low=A_min, high=A_max))                           # Random array for the Area of the Bottom Bar in the Truss
            AT1[i][j]=round(np.random.uniform(low=A_min, high=A_max))                           # Random array for the Area of the Top Bar in the Truss             
            AD11[i][j]=round(np.random.uniform(low=A_min, high=A_max))                          # Random array for the Area of the Bar going Downwards in the Truss    
            AD12[i][j]=round(np.random.uniform(low=A_min, high=A_max))                          # Random array for the Area of the Bar going Upwards in the Truss
            
            Init_population[i][j]=[VP[i][j],VL[i][j],AV1[i][j],AB1[i][j],AT1[i][j],D[i][j],AD11[i][j],AD12[i][j]]
    #
    # End
    #

    #
    # Modify the Initial population array based on the Diagonal Value
    #
    for i in range(population_size):   
        for k in range(number_of_span):
            if Init_population[i][k][5]==1.0:
                Init_population[i][k][7]=0
            elif Init_population[i][k][5]==2.0:
                Init_population[i][k][6]=0

        for l in range(3,total_number_of_parameters):
            Init_population[i][number_of_span][l]=0 
    #
    # End
    #   

    return Init_population
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

population=initial_population(max_length_of_bar,max_height_of_bar,number_of_span,A_min,A_max,population_size)   # Function call to get the Initial population

print('Initial Population \n',population)                                                                       # Printing the Initial Population         

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

nodes=np.zeros(2*(number_of_span+1))            # Initial array for nodes 
values=np.zeros((2,2*(number_of_span+1)))       # Initial array for node values 

x_coordinates=np.zeros((2*(number_of_span+1)))  # Initial array for X coordinate points 
y_coordinates=np.zeros((2*(number_of_span+1)))  # Initial array for Y coordinate points 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def node_assignments(population):
    '''
    ======================================================================
    Node dictionary assignment function
    ----------------------------------------------------------------------
    population_size: Total number of population
    ----------------------------------------------------------------------
    return: 
    node: Dictionary node and the position of the node
    n_dof: Total number of DOFs
    nodes: Array of nodes
    ======================================================================
    '''
    nodes=np.arange(2*(number_of_span+1))   # Arrange the number of nodes from 0 in the truss

    n_dof=2*len(nodes)                      # Total number of DOFs
    
    #
    # Calculating the x and y coordinates from the input array
    #
    for i in range(number_of_span+1):
        x_coordinates[2*i]=population[i][0]+max_length_of_bar*i
        x_coordinates[(2*i)+1]=x_coordinates[2*i]

        y_coordinates[2*i]=population[i][0]
        y_coordinates[(2*i)+1]=population[i][1]
    #
    # End
    #
    
    values=np.stack((x_coordinates, y_coordinates), axis=0) # Stacking the coordinate values to the values array

    node=dict(zip(nodes,zip(*values)))                      # Dictionary of the node with node number and coordinate values (x,y)

    return node,n_dof,nodes,values
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
def setup(population,node,restrained_dofs_array,force,n_dof,nodes,ele_density,ele_stiffness):
    '''
    ======================================================================
    Initial setup for FEM
    ----------------------------------------------------------------------
    population: Total number of Population
    node: Total number of nodes in each Truss 
    restrained_dofs_array: Restrained DOFs in the Truss
    force: Force vector applied in the Truss
    n_dof: Total number of DOFs
    nodes: Dictionary of node value and the coordinates of each node
    ----------------------------------------------------------------------
    return: 
    node: Total number of nodes in each Truss
    elements: Element Values in each Truss
    areas_array: Array of cross-sectional area of each element in the Truss
    degrees_of_freedom: DOFs array
    n_dof: Total number of DOFs
    densities: Dictionary of density of each element in the Truss
    stiffness: Dictionary of Stiffness of each element in the Truss
    restrained_dofs_array: Array of restrained DOFs of truss
    force_vector: Force array in the truss
    number_of_span: Total number of span in the truss
    element_number: Dictionary of element with the node route 
    ======================================================================
    '''
    print('Chromosome \n', population)              # Printing the each population 

    force_vector=dict(zip(node,force))              # Dictionary of force vector and node

    length_of_elements_array= (number_of_span*3)+1  # calculating the length of the element arrays

    #
    # Calculating the length of the element from the input array with position of 5 th element
    #
    for i in range(number_of_span+1):
        if population[i][5] == 1.0 or population[i][5]==2.0:
            length_of_elements_array+=1
        elif population[i][5]==0:
           break
        else:
            length_of_elements_array+=2
    #
    # End
    #
    element_number=np.arange(length_of_elements_array+1)    # Arranging the element number 

    node_route=[]                                           # Initial node route array

    #
    # Getting the node route from each population
    #
    for i in range((number_of_span)+1):
        if(population[i][1]!=0):
            node_route.append([nodes[2*i],nodes[(2*i)+1]])
            if(population[i][3]!=0):
                node_route.append([nodes[(2*i)],nodes[(2*i)+2]])
            if(population[i][4]!=0):
                node_route.append([nodes[(2*i)+1],nodes[(2*i)+3]])
            if(population[i][5]==1.0):
                node_route.append([nodes[(2*i)+1],nodes[(2*i)+2]])
            if(population[i][5]==2.0):
                node_route.append([nodes[(2*i)],nodes[(2*i)+3]])
            if(population[i][5]==3.0):
                node_route.append([nodes[(2*i)+1],nodes[(2*i)+2]])
                node_route.append([nodes[(2*i)],nodes[(2*i)+3]])
    #
    # End
    #

    elements=dict(zip(element_number, node_route))  # Elements dictionary with element number and route 

    areas=[]                                        # Initial array for cross-sectional areas
    #
    # Calculating the areas of the individual element from the input array
    #
    for i in range((number_of_span)+1):
        for j in [x for x in range(2,total_number_of_parameters) if x!=5]:
            if (population[i][j]!=0):
                areas.append([population[i][j]])
    #
    # End
    #
    areas_array=np.array(areas)
    
    dofs=[]                     # Initial array for DOFS
    #
    # Getting the DOfs of each node in each population
    #
    for i in range(len(node)):
        dofs.append([(2*i)+1,(2*i)+2])
    #
    # End
    #
    degrees_of_freedom=dict(zip(node,dofs)) # Dictionary of dof and node
    
    element_densities=ele_density*np.ones((len(elements)))        # Densities of each element in the truss

    densities=dict(zip(element_number,element_densities))   # Dictionary of element number and densities
    
    element_stiffness=ele_stiffness*np.ones((len(elements)))       # Stiffness of each element in the truss

    stiffness=dict(zip(element_number,element_stiffness))   # Dictionary of element number and stifness
    
    assert len(element_densities)== len(element_stiffness) == len(areas)== len(elements) # Assertion checking with size of Densities, size of Stiffness and size of elements
    assert len(restrained_dofs_array) < len(dofs)                                        # Assertion checking with size of Restrained DOFs and size of Dofs
    assert len(force_vector)== len(node)                                                 # Assertion checking with size of Force vector and size of node

    return { 'node':node,'elements':elements, 'areas':areas_array, 'degrees_of_freedom':degrees_of_freedom, \
             'ndof':n_dof,'densities' :densities, 'stiffness' :stiffness, \
             'restrained_dofs' :restrained_dofs_array,'force_vector' : force_vector, \
             'number_of_span' :number_of_span, 'element_number' :element_number }
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def node_graph(node):
    '''
    ======================================================================
    Ploting the node
    ----------------------------------------------------------------------
    node: Total number of nodes in each Truss 
    ======================================================================
    '''
    x=[i[0] for i in node.values()]             # X values from the node dictionary
    y=[i[1] for i in node.values()]             # Y values from the node dictionary
    plt.scatter(x,y,s=100,c='yellow',zorder=5)  # Scatter plot 
    offset=100/1000
    #
    # Printing the node number in the graph
    #
    for i,j in enumerate(zip(x,y)):
        plt.annotate(i+1,(j[0]-offset,j[1]-offset),zorder=10)
    #
    # End
    #
    for i_x, i_y in zip(x, y):
        plt.text(i_x, i_y, '({}, {})'.format(round(i_x,3), round(i_y,3)),va='bottom',position=(i_x+0.5,i_y+0.5))   
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
from_to_node=[]
def points(element, properties):
    '''
    ======================================================================
    Getting the Points function
    ----------------------------------------------------------------------
    element: Total elements in the Truss
    properties: function call
    ----------------------------------------------------------------------
    return: 
    from_point: Coordinate value of from node
    from_node: Value of from node of each element
    to_point: Coordinate value of to node
    to_node: Value of to node of each element
    dof: DOFs array 
    ======================================================================
    '''
    nodes=properties['node']                                # Getting nodes from Setup function call
    elements=properties['elements']                         # Getting elements from Setup function call
    degrees_of_freedom= properties['degrees_of_freedom']    # Getting degress_of_freedom from Setup function call

    from_node=elements[element][0]                          # Assigning from node to the 0 th position of the element dict
    to_node=elements[element][1]                            # Assigning to node to the 1 th position of the element dict

    from_point=np.array(nodes[from_node])                   # Assigning from point to the from node of the node dict
    to_point=np.array(nodes[to_node])                       # Assigning to point to the to node of the node dict

    
    start=degrees_of_freedom[from_node]                     # Assigning DOFS of from node
    end=degrees_of_freedom[to_node]                         # Assigning DOFS of to node
    dof=np.concatenate((start,end))

    return from_point,from_node,to_point, to_node,dof
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def draw_line(from_point, to_point, element, areas):
    '''
    ======================================================================
    Drawing line for each element 
    ----------------------------------------------------------------------
    from_point: Coordinate value of from node
    to_point: Coordinate value of to node
    element: Total elements in the Truss
    areas: Array of cross-sectional area of each element in the Truss 
    ======================================================================
    '''
    x_1=from_point[0]   # X corridinate of from point
    y_1=from_point[1]   # Y corridinate of from point
    x_2=to_point[0]     # X corridinate of to point
    y_2=to_point[1]     # Y corridinate of to point

    # Plotting line for Each X and Y coordinate of From and to point
    plt.plot([x_1,x_2],[y_1,y_2], color='green', linestyle='-', linewidth=0.5*areas[element], zorder=1)
    plt.text((x_1+x_2)/2.0, (y_1+y_2)/2.0, str(element), color="red", fontsize=12, ha='left')
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
def get_matrices(properties):
    '''
    ======================================================================
    Calculating matices
    ----------------------------------------------------------------------
    properties: Function call
    ----------------------------------------------------------------------
    return: 
    M: Global Mass matrix
    K: Global Stiffness matrix
    F: Global Force Vector
    ======================================================================
    '''
    node=properties['node']                 # Getting nodes from Setup function call
    elements=properties['elements']         # Getting elements from Setup function call
    areas=properties['areas']               # Getting areas from Setup function call
    n_dof=properties['ndof']                # Getting ndofs from Setup function call
    force_vector=properties['force_vector'] # Getting force vector from Setup function call

    K=np.zeros((n_dof,n_dof))               # Initial array for Stiffness matrix
    M=np.zeros((n_dof,n_dof))               # Initial array for Mass matrix

    node_graph(node)                        # Function call for plotting the Graph

    #
    # Caculating Elemental matrices
    #
    for element in elements:

        from_point,from_node,to_point, to_node,dof= points(element, properties)             # Getting the Node and points of each element in the truss

        draw_line(from_point, to_point, element, areas)                                     # Function call to draw elements
        x_cosine,y_cosine=dir_cosine(from_point,to_point)                                   # Function call to get directional cosine of each element
        
        length= math.sqrt((to_point[0]-from_point[0])**2+(to_point[1]-from_point[1])**2)    # Calculating length of each element

        rho=properties['densities'][element]    # Getting densities from Setup function call for each element
        area=properties['areas'][element]       # Getting cross-sectional area from Setup function call for each element
        E=properties['stiffness'][element]      # Getting stiffness from Setup function call for each element

        ck=E*area/length
        cm=rho*area*length/6.0

        m=np.array([[2,1],[1,2]])
        k=np.array([[1,-1],[-1,1]])

        tau=transformation_matrix(x_cosine,y_cosine)    # Transformation matrix function call for each element

        m_rot=tau.T.dot(m).dot(tau) 
        k_rot=tau.T.dot(k).dot(tau)

        index=dof-1
        B=np.zeros((4,n_dof))
        #
        # B matrix for calculating Rotational matrix
        #
        for i in range(4):
            B[i,index[i]]=1.0
        #
        # End
        #
        k_r_G=B.T.dot(k_rot).dot(B)
        m_r_G=B.T.dot(m_rot).dot(B)

        K+=ck*k_r_G
        M+=cm*m_r_G
    #
    # End
    #
    F=[]    # Initial array for Force array
    #
    # Getting the values of force 
    #
    for i in force_vector.values():
       F.extend(i)
    #
    # End
    #
    F=np.array(F)

    remove_index=np.array(properties['restrained_dofs'])-1  # Array to get the remove index from the restricted Dofs from Setup function call
    #
    # Reduced Mass and Stiffness matrices
    #
    for i in [0,1]:
        M=np.delete(M, remove_index, axis=i)
        K=np.delete(K, remove_index, axis=i)
    #
    # End
    #
    F=np.delete(F, remove_index)    # Reduced Force Vector

    return M,K,F
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_stresses(properties,X):
    '''
    ======================================================================
    Calculating Stress and Strains 
    ----------------------------------------------------------------------
    properties: Function call
    X: Global Displacement
    ----------------------------------------------------------------------
    return: 
    stresses_array: Stress of each element in the truss
    strains_array: Strains of each element in the truss
    lengths_array: Length of each element in the truss
    ======================================================================
    '''
    elements=properties['elements']                         # Getting elements from Setup function call
    E=properties['stiffness']                               # Getting Stiffness from Setup function call
    degrees_of_freedom=properties['degrees_of_freedom']     # Getting DOFs from Setup function call

    stresses=[] # Initial array for Stress 
    strains=[]  # Initial array for Strains 
    lengths=[]  # Initial array for Length 
    #
    # Caculating Elemental Stiffness
    #
    for element in elements:
        from_point,from_node,to_point, to_node,dof= points(element, properties)                     # Getting the Node and points of each element in the truss
        from_to_node.append([from_node,to_node])                                                    # Calculating nodal points 
        x_cosine,y_cosine=dir_cosine(from_point,to_point)                                           # Function call to get directional cosine of each element

        length= math.sqrt((to_point[0]-from_point[0])**2+(to_point[1]-from_point[1])**2)            # Calculating length of each element

        tau=transformation_matrix(x_cosine,y_cosine)                                                # Transformation matrix function call for each element
        dof_values=np.concatenate([degrees_of_freedom[from_node],degrees_of_freedom[to_node]])-1    # DOF values for from node and to node
        global_disp=[]                                                                              # Initial array for global displacement
        #
        # Global displacement calculation
        #
        for i in dof_values:
            global_disp=np.append(global_disp,X[i])
        #
        # End
        #
        q=tau.dot(global_disp)
       
        strain=(q[1]-q[0])/length   # Strain of each element
        strains.append(strain)
        
        stress=E[element]*strain    # Stress of each element
        stresses.append(stress) 
        lengths.append(length)
        stresses_array=np.array(stresses) 
        strains_array=np.array(strains)
        lengths_array=np.array(lengths)
    #
    # End
    #
    return stresses_array, strains_array,lengths_array
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def fitness_function(properties,lengths,IP,P2,stresses,strains,max_stress,max_strain,value):
    '''
    ======================================================================
    Fitness calculation for each gene in Population
    ----------------------------------------------------------------------
    properties: Function call
    lengths_array: Length of each element in the truss
    IP: Value for Genetic Algorithm penalty function
    P2: Value for Genetic Algorithm penalty function
    stresses: Stress of each element in the truss
    strains: Strains of each element in the truss
    max_stress: Allowable stress in the Truss (N/mm**2)
    max_strain: Allowable strain in the Truss (mm) 
    value: User defined values for the Constraints
    ----------------------------------------------------------------------
    return: 
    pf_array_sum: Sum of Penalty function array
    Total_weight: Total weight of each gene in the population
    Max_stress: Maximum stress of each gene
    Max_strain: Maximum strain of each gene
    ======================================================================
    '''
    elements=properties['elements']                   # Getting elements from Setup function call
    weights=np.zeros((population_size,len(elements))) # Initial Weight arrays
    pf_array_sum=np.zeros(population_size)            # Initial Penalty function array
    Max_stress=np.zeros(population_size)              # Initial maximum stress array in each population
    Max_strain=np.zeros(population_size)              # Initial maximum strain array in each population
    Total_weight=np.zeros(population_size)            # Initial total weight in each population
    #
    # Calcuating penalty function in each population based on constraints
    #
    for i in range(population_size):
        #
        # Calcuating Weight in each element
        #
        for element in elements:
            area=properties['areas'][element]
            rho=properties['densities'][element]
            weights[i][element]=(rho*area*lengths[element])
        #
        # End
        #
        #
        # Calcuating Stressa and Strain limit in each population based on constrain values
        #
        if value==1:
            maximum_stress=np.max(stresses)
            Max_stress[i]=maximum_stress
            str_lt=((maximum_stress-max_stress)/max_stress)**2 # (Sig_max-allowable_stress / allowable_Stress)**2
        elif value==2:
            maximum_strain=np.max(X)
            Max_strain[i]=max_strain
            str_lt=((maximum_strain-max_strain)/max_strain)**2 # (Str_max-allowable_strain / allowable_strain)**2
        else:
            maximum_stress=np.max(stresses)
            maximum_strain=np.max(X)
            Max_stress[i]=maximum_stress
            Max_strain[i]=maximum_strain
            str_lt=(((maximum_stress-max_stress)/max_stress)**2)+(((maximum_strain-max_strain)/max_strain)**2)
        #
        # End
        #
        pf_const=(IP+P2*(str_lt))           
        pf=[i*pf_const for i in weights]    # Calculating penalty function based on the constrain
        pf_array=np.array(pf)               
        pf_array_sum[i]=(np.sum(pf_array))  # Sum of the penalty array
        Total_weight[i]=np.sum(weights)     # Sum of the weight of each element
    return pf_array_sum,Total_weight,Max_stress,Max_strain
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
def roulette_wheel_pop(population, pf_array_sorted_rank, number):
    '''
    ======================================================================
    Parent Seletion ( Roulette wheel )
    ----------------------------------------------------------------------
    population: Total population
    pf_array_sorted_rank: Rank of pealty function with sorted
    number: Number of parent to be seleted for Genetic Algorithm
    ----------------------------------------------------------------------
    return: 
    chosen_array: Array of Selected parents
    ======================================================================
    '''
    chosen = []                         # Initial Array for selecting parent
    max_value=max(pf_array_sorted_rank)
    #
    # Randomly select parent baed on the random value
    #
    for n in range(number):
        r = np.random.uniform(0,max_value)
        #
        # Selecing parent based on the random value in the population
        #
        for (i, individual) in enumerate(population):
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
        for j in range(3):
            if(random_mutation_boolean[i][j][0]):
                population[i][j][0]=round(np.random.uniform(low=0, high=max_height_of_bar-max_length_of_bar))
            if(random_mutation_boolean[i][j][1]):
                population[i][j][1]=round(np.random.uniform(low=max_length_of_bar, high=max_height_of_bar))
            if(random_mutation_boolean[i][j][2]):
                population[i][j][2]=round(np.random.uniform(low=A_min, high=A_max))
            if(random_mutation_boolean[i][j][3]):
                population[i][j][3]=round(np.random.uniform(low=A_min, high=A_max))
            if(random_mutation_boolean[i][j][4]):
                population[i][j][4]=round(np.random.uniform(low=A_min, high=A_max))
            if(random_mutation_boolean[i][j][5]):
                population[i][j][5]=round(np.random.randint(low=1.0, high=4.0))
            if(random_mutation_boolean[i][j][6]):
                population[i][j][6]=round(np.random.uniform(low=A_min, high=A_max))
            if(random_mutation_boolean[i][j][7]):
                population[i][j][7]=round(np.random.uniform(low=A_min, high=A_max))
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
        #
        # End
        #
        X=np.array(X)

         
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
    #
    # Create a area value based on the diagonal value in the population with in the limits
    #
    for i in range(population_size):   
        for k in range(number_of_span):
            if mutated_population[i][k][5]==1.0:
                if mutated_population[i][k][6]==0:
                    mutated_population[i][k][6]=round(np.random.uniform(low=A_min, high=A_max))
                mutated_population[i][k][7]=0
            elif mutated_population[i][k][5]==2.0:
                mutated_population[i][k][6]=0
                if mutated_population[i][k][7]==0:
                    mutated_population[i][k][7]=round(np.random.uniform(low=A_min, high=A_max))
            elif mutated_population[i][k][5]==3.0:
                if mutated_population[i][k][6]==0:
                    mutated_population[i][k][6]=round(np.random.uniform(low=A_min, high=A_max))
                if mutated_population[i][k][7]==0:
                    mutated_population[i][k][7]=round(np.random.uniform(low=A_min, high=A_max))

        for l in range(3,total_number_of_parameters):
            mutated_population[i][number_of_span][l]=0 
    #
    # End
    # 
    print('Mutated Population\n', mutated_population)   # Print population after applying GA operators
    return mutated_population
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
Points=[]   # Initial array for getting the Penalty function value of each gene
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
    print('Target Value', max_stress)
    print('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
    i=0
    for x in population:
        print("Chromosome #",i, ':\n',x, '\n| Fitness: ',pf_array_sorted[i], '\n Weight: ',total_weight,'\n Maximum_stress: ',maximum_stress[j])
        print('------------------------------------------------------------')
        i=i+1
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
start = time.time()

restrained_dofs_array=get_restrained_dofs() # Getting the restrained DOFs 
force=force_vector_calc(number_of_node)     # Force vector 
pf_value_array=[]                           # Initial array for penalty function


#
# Getting the Fittness value for the first generation of the population
#
for i in range(population_size):
    node,n_dof,nodes,values=node_assignments(Init_population[i])                           # Node assignment function call
    properties=setup(population[i],node,restrained_dofs_array,force,n_dof,nodes,ele_density,ele_stiffness)    # Function call of Properties

    element=properties['elements']  # Getting elements from Setup function call
    areas=properties['areas']       # Getting cross-sectional area from Setup function call
    
    M,K,F=get_matrices(properties)  # Getting matrices function call
    X=np.linalg.inv(K).dot(F)       # Global displacement Calculation 

    print('Global displacement for chromosome '+str(i)+ ' \n',X)

    remove_index=restrained_dofs_array-1 # Restrained dofs 
    #
    # Inserting 0 in the restrained constraints in the Global displacement
    #
    for j in remove_index:
        X=np.insert(X,j,0)
    #
    # End
    #
    
    stresses, strains,lengths=get_stresses(properties,X) # Stesses and strains 

    print('Stresses for chromosome '+str(i)+'\n',stresses)
    print('strains for chromosome '+str(i)+'\n',strains)
#
# End
#
    

pf,total_weight,maximum_stress,maximum_strain=fitness_function(properties,lengths,IP,P2,stresses,strains,max_stress,max_strain,value) # Getting the fittness values and Penalized fitness based on constraints

pf_array_sorted,pf_array_sorted_rank=scaled_fitness(pf)                                                 # Sorting the rank and sorting the penalized fitness

population_1=np.zeros(np.shape(population))                                                             # Initial array for next population generation
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
print_population(population_1,0)    # Printing the population

generation_number=1
Points.append(pf_array_sorted[0])   # Append the penalty function value
#
# Loop for genetic algorithm parent selection, crossover and mutation
#
while generation_number-40:
    population=GA(population_1,pf)                          # GA algorithm function call 
    #
    # Getting the Fittness value for the each generation of the population
    #
    for i in range(population_size):
        node,n_dof,nodes,values=node_assignments(population[i])                                # Node assignment function call
        properties=setup(population[i],node,restrained_dofs_array,force,n_dof,nodes,ele_density,ele_stiffness)    # Function call of Properties

        element=properties['elements']                                                  # Getting elements from Setup function call
        areas=properties['areas']                                                       # Getting cross-sectional area from Setup function call
        
        M,K,F=get_matrices(properties)  # Getting matrices function call
        X=np.linalg.inv(K).dot(F)       # Global displacement Calculation 

        print('Global displacement \n',X)

        remove_index=restrained_dofs_array-1 # Restrained dofs 
        #
        # Inserting 0 in the restrained constraints in the Global displacement
        #
        for j in remove_index:
            X=np.insert(X,j,0) 
        #
        # End
        #
        stresses, strains,lengths=get_stresses(properties,X) # Stesses and strains 

        print('Stresses for '+str(i)+'\n',stresses)
        print('strains for '+str(i)+'\n',strains)

    

    pf,total_weight,Maximum_stress,Maximum_strain=fitness_function(properties,lengths,IP,P2,stresses,strains,max_stress,max_strain,value) # Getting the fittness values and Penalized fitness based on constraints
    
    pf_array_sorted,pf_array_sorted_rank=scaled_fitness(pf)                                                 # Sorting the rank and sorting the penalized fitness
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
    print_population(population_1,generation_number) # Printing the each population
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
with open('Genetic_algorithm.txt','w') as filename: #Writing the results in a file for plotting
    for item in Points:
        filename.write('%s\n' %item)

#-----------------------------------------------------------------------------------------------------
#Writing the best solution-----------------------------------------------------------
with open('best_solution.txt','w') as filename:
    filename.write('%s\n' %population_1[0])

for i in range(generation_number):
    plt.savefig('Truss Structure.jpg')
