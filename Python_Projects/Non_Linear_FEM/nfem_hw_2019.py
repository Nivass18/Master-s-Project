import numpy as np 
from random import random

#parameters
shape_fn=np.transpose(np.array([1/2,1/2])) #single gauss point si==1
print(shape_fn.shape)
der_shape_fn=np.transpose(np.array([-1/2,1/2])) #derivative of the shape fn
weight=2 #weight of the gauss point 
no_of_ele=3#number of elements
ri=25 #inner radius
ro=100 #outer radius
ele_len=(ro-ri)/no_of_ele #element length
poisson_ratio=0.3 
E=0.1 #young's modulus

"""
    ----------------------
    Calculating Mu and lam
    ---------------------- 
"""
mu=E/(2*(1+poisson_ratio)) #from the given problem sheet
lam=poisson_ratio*E/((1-2*poisson_ratio)*(1+poisson_ratio)) #from the given problem sheet

#Initial values of the K_global and F_ext_global as Zeros
K_global=np.zeros((no_of_ele+1,no_of_ele+1))
F_ext_global=np.zeros((no_of_ele+1,1))

#nodal positions of equal length
nodal_pstn=np.linspace(ri,ro,no_of_ele+1)

# #Material Routine
def Material_fn(lam,mu):
    C=np.array([lam+2*mu,lam,lam,lam,lam+2*mu,lam,lam,lam,lam+2*mu]).reshape(3,3)
    return C


# Element Routine
def Element_fn(nodal_pstn,der_shape_fn,shape_fn,i):

    pstn=np.array([nodal_pstn[i],nodal_pstn[i+1]])

    jacobian=np.dot(der_shape_fn,pstn)
 
    jacobian_inverse=1/jacobian

    r=np.dot(shape_fn,pstn)

    B=np.array([
                der_shape_fn*jacobian_inverse,
                [shape_fn[0]/r,shape_fn[1]/r],
                [shape_fn[0]/r,shape_fn[1]/r]
               ])
  
    B_T=B.transpose()
    print(B.shape)
    K_ele=weight*np.dot(B_T,np.dot(Material_fn(lam,mu),B))*np.power(np.dot(shape_fn,pstn),2)*jacobian
    print(np.dot(B_T,np.dot(Material_fn(lam,mu),B)))

    f_ext_ele=np.power(np.dot(np.transpose(shape_fn),pstn),2)*((2*mu*der_shape_fn*jacobian_inverse)+lam*((der_shape_fn*jacobian_inverse))+[shape_fn[0]/np.dot(np.transpose(shape_fn),pstn),shape_fn[1]/np.dot(np.transpose(shape_fn),pstn)]+[shape_fn[0]/np.dot(np.transpose(shape_fn),pstn),shape_fn[1]/np.dot(np.transpose(shape_fn),pstn)])*np.transpose(shape_fn)
    
    
    return K_ele , f_ext_ele
    # return f_ext_ele
   

for i in range(no_of_ele):
    Ae=np.zeros((2,no_of_ele+1))
    Ae[0,i]=1
    Ae[1,i+1]=1 
    AeT=Ae.transpose()

    K_ele,f_ext_ele=Element_fn(nodal_pstn,der_shape_fn,shape_fn,i)
    # print(f_ext_ele)
    K=np.dot(AeT,np.dot(K_ele,Ae)) # global K 
    F_ext=AeT@np.transpose(f_ext_ele) #global F matrix

    
   
    K_global=np.add(K_global,K)
    F_ext_global=np.add(F_ext_global,F_ext)


# # print(np.dot(F_ext_global,u))
print(K_global)
# print(F_ext_global)








