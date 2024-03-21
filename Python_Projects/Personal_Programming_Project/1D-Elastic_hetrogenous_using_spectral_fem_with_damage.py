"""
1D wave propagation in a Heterogeneous medium with force applied on the center of the structure

The values of the GLL points are obtained from the website https://en.wikipedia.org/wiki/Gaussian_quadrature and https://keisan.casio.com/exec/system/1280801905 (Lines [28 to 66])

Damage has been introduced in the structure randomly and the true length is calculated

After the damage introduction damage vector is calculated with the changes in the stiffness matrix, which has the information about the damage

The estimated value of the damage vector is calculated and the rmse error between the true and estimated damage vector is calculated, the location of the damage is calculated in which the rmse error is maximum

"""

import numpy as np 
import math 
import matplotlib.pyplot as plt
import random 

#from user

n_ele=100 #number of elements
t_len=8000 #total length (m)
rho=2000  #density of the material (kg/m3)
velocity=2500 #velocity in (m/s)
n=2 #order of polynomials
time_step = 7500 #number of time steps
dom_per=0.4 #dominant period of ricker wavelet source
alpha=0.2 #for defected element K matrix

#function to return the GLL points and the corresponding weights

def GLL_points_and_weights(n):
    #based on the n value corresponding GLL points and Weights will be returned

    if n==1:
        points = np.array([-1,1])
        weights = np.array([1,1])

    elif n==2:
        points = np.array([-1,0,1])
        weights = np.array([0.33333333, 1.33333333, 0.33333333])
    
    elif n==3:
        points = np.array([-1,-0.447213595499957939282,0.447213595499957939282,1])
        weights = np.array([0.1666666666666666666667,0.8333333333333333333333,0.833333333333333333333,0.1666666666666666666667])
    
    elif n==4:
        points =np.array([-1,-0.6546536707079771437983,0,0.6546536707079771437983,1])
        weights = np.array([0.1,0.544444444444444444444,0.7111111111111111111111,0.544444444444444444444,0.1])
    
    elif n==5:
        points = np.array([-1,-0.765055323929464692851,-0.2852315164806450963142,0.2852315164806450963142,0.765055323929464692851,1])
        weights = np.array([0.06666666666666666666667,0.3784749562978469803166,0.5548583770354863530167,0.5548583770354863530167,0.3784749562978469803166,0.06666666666666666666667])
    
    elif n==6:
        points =np.array([-1,-0.830223896278566929872,-0.4688487934707142138038,0,0.468848793470714213804,0.830223896278566929872,1])
        weights = np.array([0.04761904761904761904762,0.276826047361565948011,0.4317453812098626234179,0.487619047619047619048,0.431745381209862623418,0.2768260473615659480107,0.04761904761904761904762])
    
    elif n==7:
        points =np.array([-1,-0.8717401485096066153375,-0.5917001814331423021445,-0.2092992179024788687687,0.2092992179024788687687,
        0.5917001814331423021445,0.8717401485096066153375,1])
        weights =np.array([0.03571428571428571428571,0.210704227143506039383,0.3411226924835043647642,0.4124587946587038815671,0.412458794658703881567,0.341122692483504364764,0.210704227143506039383,0.03571428571428571428571])
    
    elif n==8:
        points =np.array([-1,-0.8997579954114601573124,-0.6771862795107377534459,-0.3631174638261781587108,0,0.3631174638261781587108,0.6771862795107377534459,0.8997579954114601573124,1])
        weights =np.array([0.02777777777777777777778,0.1654953615608055250463,0.274538712500161735281,0.3464285109730463451151,0.3715192743764172335601,0.3464285109730463451151,0.2745387125001617352807,0.165495361560805525046,0.02777777777777777777778])
    
    elif n==9:
        points =np.array([-1,-0.9195339081664588138289,-0.7387738651055050750031,-0.4779249498104444956612,-0.1652789576663870246262,0.1652789576663870246262,0.4779249498104444956612,0.7387738651055050750031,0.9195339081664588138289,1])
        weights =np.array([0.02222222222222222222222,0.1333059908510701111262,0.2248893420631264521195,0.2920426836796837578756,0.3275397611838974566565,0.3275397611838974566565,0.292042683679683757876,0.224889342063126452119,0.133305990851070111126,0.02222222222222222222222])

    return points,weights


#derivative of legrange polynomials and legendre polynomials

def lagrange(n,i,x):

    #getting the values of the GLL points and the corresponding weights
    [x_points,weights] = GLL_points_and_weights(n)
    val=1
    for j in range(-1,n):
        if j!=i:
            val=val*((x-x_points[j+1])/(x_points[i+1]-x_points[j+1]))
    return val

#inside the element with legendre polynomials with GLL points

def Legendre(n,x):
    if (n==0):
    	return x*0+1.0
    elif (n==1):
	    return x
    else:
        return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n



def lagrange_derivative(n):
    #calculation of D matrix

    #getting the values of the GLL points and the corresponding weights
    [x_points,weights] = GLL_points_and_weights(n)

    d_matrix=np.zeros([n+1,n+1]) #for the calculation of the Lagrange derivatives
    val=np.zeros([n+1,n+1]) #returns the value from Lagrange derivative

    #Values of the d matrix
    for i in range(-1,n):
        for j in range(-1,n):
            if i!=j:
                d_matrix[i+1,j+1]=Legendre(n,x_points[i+1])/Legendre(n,x_points[j+1]) * 1.0/(x_points[i+1]-x_points[j+1])
            elif i==-1 and j==-1:
                d_matrix[i+1,j+1]=-1.0/4.0*n*(n+1)
            elif i==n-1 and j==n-1:
                d_matrix[i+1,j+1]=1.0/4.0*n*(n+1)
    
    #calculation of derivative
    for i in range(-1,n):
        for j in range(-1,n):
            sum=0
            for k in range(-1,n):
                sum=sum+d_matrix[j+1,k+1]*lagrange(n,i,x_points[k+1])
            val[i+1,j+1]=sum
        
    return(val)


lagrange_1_derivative=lagrange_derivative(n) #1st derivative of the Lagrange polynomial is assigned

#Force vector calculation

def gaussian_force(global_dt,dom_per):
    a=int(2*dom_per/global_dt)
    res=np.zeros(a)
    t_0=dom_per/global_dt
    value=4/dom_per
    for i in range(0,a):
        t=((i+1)-t_0)*global_dt
        res[i]=-2*value*t*math.exp(-(value*t)**2)
    return res

#Element level matrices function

#function to calculate the element mass matrix

def ele_mass_matrix_fn(weights,rho,jacobian):
    for i in range(0,n+1):
        ele_mass_matrix[i] = weights[i]*rho[i]*jacobian # mass matrix diagonal element is stored as an 1d array
    return ele_mass_matrix

#Function to calculate the element stiffness matrix

def ele_stiffness_matrix_fn(mu,inv_jacobian,weights,l):
    for i in range (-1,n):
            for j in range(-1,n):
                sum=0
                for k in range(-1,n):
                    sum=sum+mu[k+1+l]*weights[k+1]*inv_jacobian**2*jacobian*lagrange_1_derivative[i+1,k+1]*lagrange_1_derivative[j+1,k+1]
                ele_stiffness_matrix[i+1,j+1]=sum
    return ele_stiffness_matrix



#calculation

e_len=t_len/n_ele #element length (same length for each element)

n_global = n_ele*n+1 #global size of the stiffness and the mass matrix

random_ele_number=random.randint(0,n_ele) #random number for the defected ele

length_of_defected_ele=random_ele_number*e_len #to calculate the length of the randomly generated defected element

#jacobian and inverse jacobian calculation

jacobian=e_len/2
inv_jacobian=1/jacobian

#global space domain calculation

x_global=np.zeros(n_global)
k=(t_len/(n_global-1))
val=0
for i in range(0,n_global):
    x_global[i]=val
    val=val+k 

x_global_min=min(np.diff(x_global)) #minimum value difference of global space domain
courant_value=0.1 #courant value 
global_dt=courant_value*x_global_min/velocity #Global time step

#for hetrogeneous model with fluction in wave velocity and density at certain elements

no_of_ele=25 #number of elements in the velocity reduction
reduction=0.4 #value of the reduction of velocity

velocity=velocity*np.ones((n_ele*n+1)) #calculates the velocity array
rho=rho*np.ones((n_ele*n+1)) #calculates the density array
velocity[n_ele-no_of_ele:n_ele+no_of_ele]=max(velocity)*reduction #calculates the reduction of velocity in certain region

mu=rho*velocity**2 #shear modulus [shear velocity = sqrt(shear_modulus/density)]

#getting the values of the GLL points and the corresponding weights

[x_points,weights] = GLL_points_and_weights(n)

#initiation of matrices

ele_mass_matrix=np.zeros((n+1,1),dtype=float) #element mass matrix
mass_matrix = np.zeros(n_global,dtype=float) # mass matrix for diagonal elements calculation
global_mass_matrix=np.zeros((n_global,n_global),dtype=float) # global mass matrix
global_inv_mass_matrix=np.zeros((n_global,n_global),dtype=float) #global inverse mass matrix
ele_stiffness_matrix=np.zeros((n+1,n+1),dtype=float) #element stiffness matrix
global_stiffness_matrix=np.zeros((n_global,n_global),dtype=float) #global stiffness matrix
global_stiffness_matrix_diagonal=np.zeros((n_global,1),dtype=float) #global stiffness matrix

#Displacement, force and acceleration vector assignment with consideration of damage

u=np.zeros(n_global)
u_old=u
u_new=u
force=u
acceleration=u
#velocity=u

ele_mass_matrix=ele_mass_matrix_fn(weights,rho,jacobian)  # element mass matrix function call

#global mass matrix

#Diagonal element of the mass matrix is calculated
k=-1
for i in range(1,n_ele+1):
    for j in range(0,n+1):
        k=k+1
        if i>1:
            if j==0:
                k=k-1
        mass_matrix[k]=mass_matrix[k]+ele_mass_matrix[j]

# Arrangement of the diagonals of the mass matrix to the global mass matrix

for i in range(n_global):
    for j in range(n_global):
        #Main diagonal elements of the global mass matrix is assigned to the values
        if i==j:
            global_mass_matrix[i,j]=mass_matrix[j]
        else:
            global_mass_matrix[i,j]=0

# Inverse of a global mass matrix

for i in range(n_global):
    for j in range(n_global):
        if i==j:
            global_inv_mass_matrix[i,j]=1.0/global_mass_matrix[i,j] 
        else:
            global_inv_mass_matrix[i,j]=0


#global_striffness matrix with consideration of damage
l=0

for p in range(1,n_ele+1):
    i0=(p-1)*n+1
    j0=i0
    #Based on random ele number for damage the element stiffness matrix is selected
    if random_ele_number==0:
        ele_stiffness_matrix=ele_stiffness_matrix_fn(mu,inv_jacobian,weights,l)
    elif p==random_ele_number: #damage location is not zero
        ele_stiffness_matrix=alpha*ele_stiffness_matrix_fn(mu,inv_jacobian,weights,l) #alpha*element stiffness matrix 
    else:
       ele_stiffness_matrix=ele_stiffness_matrix_fn(mu,inv_jacobian,weights,l)
    l+=n
    for i in range(-1,n):
        for j in range(-1,n):
            global_stiffness_matrix[i0+i,j0+j]+=ele_stiffness_matrix[i+1,j+1] #assiging the element stiffness matrix to the global stiffness matrix

#source and the location of the source

source=gaussian_force(global_dt,dom_per)

source_location=int(np.floor(n_global/2)) #Location of the source

#displacement calculation with damage 

for i in range(time_step):
    force=np.zeros(n_global) #Array to store the value of the force

    #Force Calculation at the center of structure
    if i < len(source):
        force[source_location-1]=source[i-1]
    
    #Rigid boundary condition at both ends    
    #for j in range(2,(n_global-1)):
     #   u_new=global_dt**2*global_inv_mass_matrix@(force-global_stiffness_matrix@u)+2*u-u_old
    
    #Free boundary condition 
    u_new=global_dt**2*global_inv_mass_matrix@(force-global_stiffness_matrix@u)+2*u-u_old #Newmark scheme
    acceleration=(u_new-2*u+u_old)/(global_dt**2) #acceleration based on the newmark scheme
    #velocity=(u_new-u)/global_dt
    u_old,u=u,u_new



#np.savetxt("Force_at_end_Displacement_values_with_order_9.txt", u_new, fmt="%s")

#Plot for the damaged displacement
plt.plot(x_global,u_new)
plt.title("Displacement graph with consideration of damage")
plt.xlabel(" Global coordinates")
plt.ylabel("u_new")
plt.grid()
plt.show()

#plt.savefig('Force_at_end_Displacement_values_with_order_9.jpg')


#displacement arrays without damage

global_stiffness_matrix_undamaged=np.zeros((n_global,n_global),dtype=float) #global stiffness matrix without damage 


#global stifness matrix calculation without damage
l=0
for p in range(1,n_ele+1):
    i0=(p-1)*n+1
    j0=i0
    ele_stiffness_matrix=ele_stiffness_matrix_fn(mu,inv_jacobian,weights,l) #Element Stiffness matrix function call
    l+=n
    for i in range(-1,n):
        for j in range(-1,n):
            global_stiffness_matrix_undamaged[i0+i,j0+j]+=ele_stiffness_matrix[i+1,j+1] #assigning the values of element stiffness matrix to the global stiffness matrix
 
#Displacement and acceleration array initalization to store the values without consideration of damage        
u_undamaged=np.zeros(n_global)
u_new_undamaged=u_undamaged
u_old_undamaged=u_undamaged
acceleration_undamaged=u_undamaged

#displacement  calculation without damage  
  
for i in range(time_step):
    force=np.zeros(n_global) #Array to store the value of the force

    #Force Calculation at the center of structure
    if i < len(source):
        force[source_location-1]=source[i-1]

    #Rigid boundary condition at both ends   
    #for j in range(2,(n_global-1)):
        #u_new_undamaged=global_dt**2*global_inv_mass_matrix@force+2*u-u_old

    #Free boundary condition    
    u_new_undamaged=global_dt**2*global_inv_mass_matrix@(force-global_stiffness_matrix_undamaged@u_undamaged)+2*u_undamaged-u_old_undamaged #Newmark scheme
    acceleration_undamaged=(u_new_undamaged-2*u_undamaged+u_old_undamaged)/(global_dt**2) #acceleration based on the newmark scheme
    u_old_undamaged,u_undamaged=u_undamaged,u_new_undamaged 

#Printing the values
print("The number of elements ",n_ele)
print("Total length of the structure ",t_len)
print("Order of the polynomial ",n)
print("The randomly generated element of the damage",random_ele_number)
print("The element length",e_len)
print("The location of the damage",length_of_defected_ele)



#Inverse problem of finding rmse error between damage vector and location of damage prediction

#Changes in the stiffness matrix
global_stiffness_matrix_changes=np.zeros((n_global,n_global),dtype=float) #changes in global stiffness matrix 

damage_vector=np.zeros((n_global,1),dtype=float) #Damage vector which changes based on the stiffness vector change

global_stiffness_matrix_changes=global_stiffness_matrix_undamaged-global_stiffness_matrix #finds the changes in the global stifness matrix

damage_vector=global_stiffness_matrix_changes@u_new #True value of Damage Vector calculation

#Calculation of the estimated damage vector
for i in range(n_global):
    b=global_mass_matrix@(acceleration_undamaged-acceleration)+global_stiffness_matrix@(u_new_undamaged-u_new) 
    b_reshape=b.reshape(n_global,1)


difference= b-damage_vector #Difference between the Estimated and the True value of the damage vector is calculated and used for plot

rmse=np.zeros(n_global,dtype="float") #Array to store the rmse error

#calculation of the rmse error
for i in range(n_global):
    rmse[i]=np.sqrt((b[i]-damage_vector[i])**2/(n_global))

#based on Random Element number location of the damage is predicted    
if random_ele_number==0: 
    print("There is no damage in the structure and the Structure is healthy")
else:
    a=(random_ele_number-0.5)*(t_len/n_ele) #lowest limit
    b=(random_ele_number+0.5)*(t_len/n_ele) #highest limit 

    #fiding the index of the upper and lower limits
    for i in range (n_global):
        if x_global[i]==a:
            a_index=i+1
        elif a==0:
            a_index=0


    for j in range (n_global):
        if x_global[j]==b:
            b_index=j+1
        elif b>=x_global[n_global-1]:
            b_index=n_global-1


 
    #Splits the Rmse Error based on the upper and lower limit
    rmse_split=rmse[a_index:b_index+1]
    print("The changes in the RMSE Error ",rmse_split)

    maximum_value=max(rmse_split) #maximum value of the rmse error is calculated

    #index of the maximum value in the rmse found and the corresponding index value in the global coordinate is printed as location of the damage
    for k in range (b_index+1):
        if rmse[k]==maximum_value:
            print("The position of the damage",x_global[k])


#Plots for the True value of the damage vector, estimated value of the damage vector and the difference between them

plt.plot(x_global,rmse,color="orange",label="RMSE")
plt.plot(x_global,damage_vector,color="blue",label="True Damage vector")
plt.plot(x_global,b_reshape,color="black",label="Estimated Damage vector")
plt.title("True and estimated Value between damage vector and rmse error")
plt.xlabel("Global coordinates")
plt.ylabel("Ture and estimated value of damage vector, rmse")
plt.legend(loc="best")
plt.grid()
plt.show()
