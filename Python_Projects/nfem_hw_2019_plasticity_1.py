import matplotlib.pyplot as plt
import numpy as np

#Time step parameters
time_step=0.01
time_step_array=np.arange(time_step,1,time_step)

#INPUT PARAMETERS   
poissons_ratio=0.3
E=0.1
volumetric_strain=-0.01
yield_stress=100e-6

#Mu and Lamda
mu=E/(2*(1+poissons_ratio))
Lamda=poissons_ratio*E/((1-2*poissons_ratio)*(1+poissons_ratio))

#MESH PARAMETERS
inner_radius = 25
outer_radius = 100
no_of_ele = 10
meshrefinementfactor = 5 #ratio of element sizes at outer and inner radius

#MESH REFINEMENT
q = meshrefinementfactor**(1/(no_of_ele-1))
l_ele = (outer_radius-inner_radius)*(1-q)/(1-meshrefinementfactor*q)
r_node = inner_radius
coordinate=np.array([inner_radius])

#MESH GENERATION
for i in range(no_of_ele):
    r_node = r_node +l_ele
    coordinate = np.append(coordinate,r_node)
    l_ele=l_ele*q

C = np.array([[Lamda+2*mu,Lamda,Lamda],
              [Lamda,Lamda+2*mu,Lamda],
              [Lamda,Lamda,Lamda+2*mu]])
        
#Initial values of the K_global and F_ext_global as Zeros
K_global=np.zeros((no_of_ele+1,no_of_ele+1))
Global_epsilon_p=np.zeros((no_of_ele,3,1))
F_int_global=np.zeros((no_of_ele+1,1))

#Initial Guess of Displacement for the Newton Raphson Method
u=np.linspace(0,(1/3*time_step_array[0]*(-volumetric_strain)*coordinate[0]),(no_of_ele+1)).reshape(no_of_ele+1,1)
global_u=np.flip(u,0)
u_reduced=np.delete(global_u,(0),axis=0)

#Plasticity Parameters
k=3*Lamda+2*mu
delta_ij=np.eye(3)
delta_ij_delta_kl=np.outer(delta_ij,delta_ij)
delta_ik_delta_jl=np.eye(9)
delta_il_delta_jk=np.zeros((9,9))
delta_il_delta_jk[0,0]=1;delta_il_delta_jk[1,3]=1;delta_il_delta_jk[2,6]=1;delta_il_delta_jk[3,1]=1;delta_il_delta_jk[4,4]=1;delta_il_delta_jk[5,7]=1;delta_il_delta_jk[6,2]=1;delta_il_delta_jk[7,5]=1;delta_il_delta_jk[8,8]=1
I_ijkl=(1/2*(delta_ik_delta_jl+delta_il_delta_jk))-(1/3*delta_ij_delta_kl)

#For plot
sigma_rr_ri=[]

# Element Routine
def Element_fn(coordinate,Lambda,mu,tau):
    der_shape_fn=np.array([-1/2,1/2])
 
    jacobian=der_shape_fn@np.array([[coordinate[i]],[coordinate[i+1]]])
    
    jacobian_inverse=np.asscalar(1/jacobian)

    B=np.array([
                [-1/2*jacobian_inverse,1/2*jacobian_inverse],
                [1/(coordinate[i]+coordinate[i+1]),1/(coordinate[i]+coordinate[i+1])],
                [1/(coordinate[i]+coordinate[i+1]),1/(coordinate[i]+coordinate[i+1])]
                ])
  
    B_T=np.transpose(B)
    
    C,new_sigma=Material_fn(coordinate,global_u)

    K_ele=2*B_T@C@B*jacobian*((coordinate[i]+coordinate[i+1])/2)**2

    f_int_ele=2*(B_T@new_sigma)*jacobian*((coordinate[i]+coordinate[i+1])/2)**2
   
    return K_ele,f_int_ele


# #Material Routine
def Material_fn(coordinate,global_u):
    der_shape_fn=np.array([-1/2,1/2])
 
    jacobian=der_shape_fn@np.array([[coordinate[i]],[coordinate[i+1]]])
    
    jacobian_inverse=np.asscalar(1/jacobian)

    B=np.array([
                [-1/2*jacobian_inverse,1/2*jacobian_inverse],
                [1/(coordinate[i]+coordinate[i+1]),1/(coordinate[i]+coordinate[i+1])],
                [1/(coordinate[i]+coordinate[i+1]),1/(coordinate[i]+coordinate[i+1])]
                ])
    
    epsilon = B @ np.array([global_u[i],global_u[i+1]])
    
    sigma_trail=C@(epsilon-Global_epsilon_p[i])

    sigma_trail_dvt = sigma_trail-(1/3*np.sum(sigma_trail))

    sigma_trail_eql=np.sqrt(3/2*(np.sum(np.square(sigma_trail_dvt))))

    sigma_trail_dvt_tensor=np.diagflat(sigma_trail_dvt)

    sigma_trail_dvt_outer=np.outer(sigma_trail_dvt_tensor,sigma_trail_dvt_tensor)

    if (sigma_trail_eql-yield_stress) < 0:
        sigma_rr[i]=np.asscalar(sigma_trail[0])
        sigma_phi[i]=np.asscalar(sigma_trail[1])
        return C,sigma_trail
    else:
        delta_lamda=(sigma_trail_eql-yield_stress)/(3*mu)

        epsilon_p[i]=Global_epsilon_p[i]+delta_lamda*1.5*(sigma_trail_dvt/sigma_trail_eql)

        new_trial_stress=C@(epsilon-epsilon_p[i])

        sigma_rr[i]=np.asscalar(new_trial_stress[0])
        sigma_phi[i]=np.asscalar(new_trial_stress[1])

        new_trial_stress_dev=new_trial_stress-(1/3*np.sum(new_trial_stress))

        new_stress_eql=np.sqrt(3/2*(np.sum(np.square(new_trial_stress_dev))))

        C_plastic = k/3*(delta_ij_delta_kl)+2*mu*((sigma_trail_eql-3*mu*delta_lamda)/sigma_trail_eql)*I_ijkl-(3*mu*(1/sigma_trail_eql**2)*sigma_trail_dvt_outer)

        value = C_plastic[0::4]
        C_t = value[np.nonzero(value)].reshape(3,3)

        return C_t,new_trial_stress

        
    

for time,tau in  enumerate(time_step_array):
    sigma_rr=np.zeros_like(coordinate)
    sigma_phi=np.zeros_like(coordinate)
    global_u[0]=1/3*tau*(-volumetric_strain)*inner_radius
    delta_u=np.array([1,1])
    G_matrix_red = np.array([1,1])
    epsilon_p=np.zeros((no_of_ele,3,1))
    while np.linalg.norm(delta_u,np.inf)>(0.005*np.linalg.norm(global_u,np.inf)) or np.linalg.norm(G_matrix_reduced,np.inf)>(0.005*np.linalg.norm(F_int_global,np.inf)):
        K_global=np.zeros((no_of_ele+1,no_of_ele+1))
        F_int_global=np.zeros((no_of_ele+1,1))
        for i in range(no_of_ele):
            K_ele,f_int_ele = Element_fn(coordinate,Lamda,mu,tau)
            Ae=np.zeros((2,no_of_ele+1))
            Ae[0,i]=1
            Ae[1,i+1]=1 
            AeT=np.transpose(Ae)
                    
            
            #Global K Matrix         
            K=AeT@K_ele@Ae  
            K_global=np.add(K_global,K)

            #Global F_int matrix
            F_int=AeT@f_int_ele 
            F_int_global=np.add(F_int_global,F_int)
            
        K_global_red=np.delete(K_global,0,axis=0)
        K_global_red=np.delete(K_global_red,0,axis=1)
       
        F_ext_global=np.zeros((no_of_ele+1,1))

        #Newton Raphson Method
        G_matrix = F_int_global-F_ext_global
        
        G_matrix_reduced = np.delete(G_matrix,(0),axis=0)
        
        delta_u=np.linalg.inv(K_global_red)@G_matrix_reduced
        u_reduced=u_reduced-delta_u
    
        global_u = np.insert(u_reduced,(0),1/3*tau*(-volumetric_strain)*inner_radius).reshape(no_of_ele+1,1)
    sigma_rr_ri.append(sigma_rr[0])
    Global_epsilon_p=epsilon_p

    
#Plots
fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(10,5))

ax[0,0].plot(coordinate,global_u,'k')
ax[0,0].set_title('Global Displacement and Coordinate')
ax[0,0].set_xlabel('Coordinates')
ax[0,0].set_ylabel('Global_displacement')
ax[0,0].grid(True)
ax[0,1].plot(coordinate,sigma_rr,'k')
ax[0,1].set_title('Sigma_rr and Coordinate')
ax[0,1].set_xlabel('Coordinates')
ax[0,1].set_ylabel('Sigma_rr')
ax[0,1].grid(True)
ax[1,0].plot(coordinate,sigma_phi,'k')
ax[1,0].set_title('Sigma_phi and Coordinate')
ax[1,0].set_xlabel('Coordinates')
ax[1,0].set_ylabel('Sigma_phi')
ax[1,0].grid(True)
ax[1,1].plot(time_step_array,sigma_rr_ri,'k')
ax[1,1].set_title('Time_step and Sigma_rr_ri')
ax[1,1].set_xlabel('Time-step')
ax[1,1].set_ylabel('Sigma_rr_ri')
ax[1,1].grid(True)
plt.tight_layout()
plt.savefig('plots.png')