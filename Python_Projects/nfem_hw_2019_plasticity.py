import numpy as np 
import matplotlib.pyplot as plt


#parameters
'''
---------------------------------------
Initial Parameter for the given problem
---------------------------------------
'''

shape_fn=np.transpose(np.array([[1/2,1/2]])) #single gauss point si==1
der_shape_fn=np.transpose(np.array([[-1/2,1/2]])) #derivative of the shape fn
weight=2 #weight of the gauss point 
no_of_ele= 2 #number of elements
ri=5 #inner radius
ro=20 #outer radius
poisson_ratio=0.2 #Poission Ratio
E=0.2 #young's modulus
epsilon_v= -0.005 #Volumetric strain
tau=np.linspace(0.1,1,10)

sigma_y = 0.0002 #yield stress 
epsilon_p=0 #Plastic Strain

delta_u=1


"""
    ----------------------
    Calculating Mu and lam
    ---------------------- 
"""
mu=E/(2*(1+poisson_ratio)) #from the given problem sheet
Lambda=((poisson_ratio*E)/((1-2*poisson_ratio)*(1+poisson_ratio))) #from the given problem sheet

#Initial values of the K_global and F_ext_global as Zeros
K_global=np.zeros((no_of_ele+1,no_of_ele+1))
F_ext_global=np.zeros((no_of_ele+1,1))

#nodal positions of equal length
nodal_pstn=np.linspace(ri,ro,no_of_ele+1).reshape(no_of_ele+1,1)

#Initial Guess of Displacement for the Newton Raphson Method
u=np.linspace(0,(1/3*(epsilon_v)*ri),(no_of_ele+1)).reshape(no_of_ele+1,1)
global_u=np.flip(u,0)
# print(global_u.shape)
u_reduced=np.delete(global_u,(0),axis=0)

#Exact Elastic Displacement 
# u_elastic=((ri)**3*(-epsilon_v)*tau)/(3*(nodal_pstn)**2)

# #Material Routine
def Material_fn(Lambda,mu,i,B):
    C = np.array([[Lambda+2*mu,Lambda,Lambda],
              [Lambda,Lambda+2*mu,Lambda],
              [Lambda,Lambda,Lambda+2*mu]])

    epsilon = B @ np.array([global_u[i],global_u[i+1]])
    sigma_trail=C@(epsilon-epsilon_p)
    Identity=np.array([1,1,1]).reshape(3,1)
    sigma_dvt = sigma_trail-(1/3*(np.trace(sigma_trail)*Identity))
    sigma_dvt_square=np.square(sigma_dvt)
    sigma_eql=np.sqrt(3/2*np.sum(sigma_dvt_square,axis=0))
    if (sigma_eql-sigma_y) < 0:
        print('Elastic ',i)
    else:
        print('Plastic ',i)
    return C
    

# Element Routine
def Element_fn(nodal_pstn,der_shape_fn,shape_fn,i):

    pstn=np.array([[nodal_pstn[i]],[nodal_pstn[i+1]]]).reshape(2,1)
 
    jacobian=np.dot(np.transpose(der_shape_fn),pstn)
    
    jacobian_inverse=np.asscalar(1/jacobian)
    
    r=(np.transpose(shape_fn)@pstn)
 
    B_1=np.matrix(der_shape_fn*jacobian_inverse)

    B=np.array([
                     [B_1[0,0],B_1[1,0]],
                    [shape_fn[0,0]/np.asscalar(r),shape_fn[1,0]/np.asscalar(r)],
                    [shape_fn[0,0]/np.asscalar(r),shape_fn[1,0]/np.asscalar(r)]
                ])
    # print(B.shape)
    B_T=B.transpose()
    
    K_ele=weight*np.dot(B_T,np.dot(Material_fn(Lambda,mu,i,B),B))*np.power(np.dot(np.transpose(shape_fn),pstn),2)*jacobian
    
    sigma_rr=(2*mu*jacobian_inverse)*((-global_u[i]+global_u[i+1])/2)+Lambda*(0.001*epsilon_v)
    # print(sigma_rr)
    f_ext_ele=np.array([-np.asscalar(sigma_rr)*nodal_pstn[i]**2,
                           np.asscalar(sigma_rr)*nodal_pstn[i+1]**2])#.reshape(2,1)
    
    # print(f_ext_ele.shape)
    return K_ele , B, f_ext_ele
    

  
for i in tau:
    global_u[0]=1/3*i*(epsilon_v)*ri
    # print(i)
    while True:
        for i in range(no_of_ele):
            Ae=np.zeros((2,no_of_ele+1))
            Ae[0,i]=1
            Ae[1,i+1]=1 
            AeT=Ae.transpose()
                    
            K_ele , B, f_ext_ele = Element_fn(nodal_pstn,der_shape_fn,shape_fn,i)
                    
            K=np.dot(AeT,np.dot(K_ele,Ae)) # global K 
            F_ext=np.dot(AeT,f_ext_ele) #global F matrix
                    
            K_global=np.add(K_global,K)
            F_ext_global=np.add(F_ext_global,F_ext)
                    # F_ext_global[:,no_of_ele]=0
        K_global_red=np.delete(K_global,0,axis=0)
        K_global_red=np.delete(K_global_red,0,axis=1)
# print(F_ext_global)       
       
        #Newton Raphson Method
        G_matrix = K_global@global_u
        
        G_matrix_reduced = np.delete(G_matrix,(0),axis=0)
        
        delta_u=np.linalg.inv(K_global_red)@G_matrix_reduced
        u_reduced=u_reduced-delta_u
         # print(u_reduced)
        global_u = np.insert(u_reduced,(0),1/3*i*(-epsilon_v)*ri).reshape(no_of_ele+1,1)
        if np.linalg.norm(delta_u) > (0.005*np.linalg.norm(u_reduced)):
            break
        # u_reduced=np.delete(global_u,(0),axis=0)
        # # initial_nodal_pstn=nodal_pstn
        # # nodal_pstn=initial_nodal_pstn+global_u
        
        # epsilon = B@(Ae@global_u)
        # print(epsilon)
        
        # # print(sigma_dvt)
        
        # print(sigma_eql)
        # Elastic_Plastic(sigma_eql,sigma_y)
        # # sigma_trail = E*(epsilon - epsilon_p)
        # # print(global_u.shape)
        # # print(global_u)
        # # print(global_u.shape)
        # # u=next_u_updated
        # #
        # # 