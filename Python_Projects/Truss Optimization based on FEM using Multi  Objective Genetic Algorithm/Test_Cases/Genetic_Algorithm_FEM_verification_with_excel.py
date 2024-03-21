'''
Basic FEM problem with 
youngs_modulus =70000 Pa; 
density=1 Kg/mm^3; 
max_length_of_bar = 3000 mm
max_height_of_bar = 3000 mm
Force [3] = -50000
force [7] = -100000
force [12] = -50000
Prescribed DOF = [1 2 10]

'''

import numpy as np
import math 
import matplotlib.pyplot as plt
import time
import sys

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
    total_number_of_parameters=int(input("Enter the number of parameters in each bar "))    # Total number of parameters in each bar
    number_of_span=int(input("Enter the number of span in the truss structure "))           # Total number of spans in the truss  

    max_length_of_bar=int(input("Enter the length of the bar "))                            # Maximum length of the bar 
    max_height_of_bar=int(input("Enter the height of the bar "))                            # Maximum height of the bar

except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
try:
    ele_density=float(input("Enter the density of element "))                                 # Density of element
    ele_stiffness=float(input("Enter the stiffness of element "))                             # Stiffness of element

except ValueError:
    print("\n Only integer values are valid")
    sys.exit(0)  

number_of_node=2*(number_of_span+1)                                                     # Total number of nodes in the Truss
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

input_array=[[0,10,8,9,10,2,0,5],[1,8,5,7,6,3,8,9],[0,9,10,0,0,0,0,0]] # For the specific Problem

print('Input array is \n ', input_array)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def get_all_elements_in_input_array(input_array):
    '''
    ======================================================================
    Get all the elements from the input array Function
    ----------------------------------------------------------------------
    input_array: Input array for the truss
    ----------------------------------------------------------------------
    return: 
    count: Total number of Elements in the truss
    ======================================================================
    '''
    count=0                     # Initial value for get the count of elements
    #
    # Counting the number of elements in the population
    #
    for element in input_array:
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
def node_assignments(input_array):
    '''
    ======================================================================
    Node dictionary assignment function
    ----------------------------------------------------------------------
    input_array: Input array for the truss
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
        x_coordinates[2*i]=input_array[i][0]+max_length_of_bar*i
        x_coordinates[(2*i)+1]=x_coordinates[2*i]

        y_coordinates[2*i]=input_array[i][0]
        y_coordinates[(2*i)+1]=input_array[i][1]
    #
    # End
    #
    
    values=np.stack((x_coordinates, y_coordinates), axis=0) # Stacking the coordinate values to the values array

    node=dict(zip(nodes,zip(*values)))                      # Dictionary of the node with node number and coordinate values (x,y)

    return node,n_dof,nodes, values
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
def setup(input_array,node,restrained_dofs_array,force,n_dof,nodes,ele_density,ele_stiffness):
    '''
    ======================================================================
    Initial setup for FEM
    ----------------------------------------------------------------------
    input_array: Input array for the truss
    node: Total number of nodes in each Truss 
    restrained_dofs_array: Restrained DOFs in the Truss
    force: Force vector applied in the Truss
    n_dof: Total number of DOFs
    nodes: Dictionary of node value and the coordinates of each node
    ele_density: Density of element
    ele_stiffness: Stiffness of element
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
    force_vector=dict(zip(node,force))              # Dictionary of force vector and node

    length_of_elements_array= (number_of_span*3)+1  # calculating the length of the element arrays

    #
    # Calculating the length of the element from the input array with position of 5 th element
    #
    for i in range(number_of_span+1):
        if input_array[i][5] == 1 or input_array[i][5]==2:
            length_of_elements_array+=1
        elif input_array[i][5]==0:
           break
        else:
            length_of_elements_array+=2
    #
    # End
    #
    element_number=np.arange(length_of_elements_array+1) # Arranging the element number

    node_route=[]                                        # Initial node route array

    #
    # Getting the node route from each population
    #
    for i in range((number_of_span)+1):
        if(input_array[i][1]!=0):
            node_route.append([nodes[2*i],nodes[(2*i)+1]])
            if(input_array[i][3]!=0):
                node_route.append([nodes[(2*i)],nodes[(2*i)+2]])
            if(input_array[i][4]!=0):
                node_route.append([nodes[(2*i)+1],nodes[(2*i)+3]])
            if(input_array[i][5]==1):
                node_route.append([nodes[(2*i)+1],nodes[(2*i)+2]])
            if(input_array[i][5]==2):
                node_route.append([nodes[(2*i)],nodes[(2*i)+3]])
            if(input_array[i][5]==3):
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
            if (input_array[i][j]!=0):
                areas.append([input_array[i][j]])
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

    element_densities=ele_density*np.ones((len(elements)))            # Densities of each element in the truss

    densities=dict(zip(element_number,element_densities))   # Dictionary of element number and densities

    element_stiffness=ele_stiffness*np.ones((len(elements)))        # Stiffness of each element in the truss

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
        plt.text(i_x, i_y, '({}, {})'.format(i_x, i_y),va='bottom',position=(i_x+0.5,i_y+0.5))
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
    plt.plot([x_1,x_2],[y_1,y_2], color='green', linestyle='-', linewidth=0.5, zorder=1)
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
        
        x_cosine,y_cosine=dir_cosine(from_point,to_point)                                   # Function call to get directional cosine of each elemen
        
        length= math.sqrt((to_point[0]-from_point[0])**2+(to_point[1]-from_point[1])**2)    # Calculating length of each element

        rho=properties['densities'][element]    # Getting densities from Setup function call for each element
        area=properties['areas'][element]       # Getting cross-sectional area from Setup function call for each element
        E=properties['stiffness'][element]      # Getting stiffness from Setup function call for each element

        ck=E*area/length # Calculating Ck value 
        cm=rho*area*length/6.0 # Calculating Cm value

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
    
    print('Determinat of K before reducing',np.linalg.det(K)) # For checking singularity before reducing
    #
    # Reduced Mass and Stiffness matrices
    #
    for i in [0,1]:
        M=np.delete(M, remove_index, axis=i)
        K=np.delete(K, remove_index, axis=i)
    #
    # End
    #
    print('Determinant of K after reducing ', np.linalg.det(K)) #for checking non-singularity after reducing

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
    
        draw_line(from_point, to_point, element, stress)
    #
    # End
    #   

    return stresses_array, strains_array,lengths_array
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
def print_results(X,stresses,strains):
    '''
    ======================================================================
    Printing Results
    ----------------------------------------------------------------------
    X: Global Displacement
    stresses: Stress in each element
    strains: Strains in each element
    ----------------------------------------------------------------------
    ======================================================================
    '''
    print("Global_displacements are \n", X)
    print("Global stresses are \n", stresses)
    print('Global strains are \n', strains)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------#----------------------------------------------------------------------------------------------------------------------------------------------------------------
start = time.time()
restrained_dofs_array=get_restrained_dofs() # Getting the restrained DOFs 
force=force_vector_calc(number_of_node)     # Force vector 
node,n_dof,nodes,values=node_assignments(input_array)                           # Node assignment function call
properties=setup(input_array,node,restrained_dofs_array,force,n_dof,nodes,ele_density,ele_stiffness)             # Function call of Properties

element=properties['elements']  # Getting elements from Setup function call
node=properties['node']         # Getting node from Setup function call
remove_index=np.array(properties['restrained_dofs'])-1 # Restrained dofs 

M,K,F=get_matrices(properties)  # Getting matrices function call

X=np.linalg.inv(K).dot(F)       # Global displacement Calculation  

#
# Inserting 0 in the restrained constraints in the Global displacement
#
for i in remove_index:
    X=np.insert(X,i,0) 
#
# End
#
    
stresses, strains,lengths=get_stresses(properties,X) # Stesses and strains 

print_results(X,stresses,strains)                   # Printing the values 

end = time.time()
print('Total time elasped',str(end - start)+str(' in seconds'))
#-----------------------------------------------------------------------------------------------------
#Writing the results in a file for plotting-----------------------------------------------------------
with open('Genetic_Algorithm_FEM_with_excel.txt','w') as filename:
        filename.write('Stress are %s\n' %stresses)
        filename.write('Global displacements are %s\n' %X)

plt.savefig('Truss Structure_with_excel.jpg')