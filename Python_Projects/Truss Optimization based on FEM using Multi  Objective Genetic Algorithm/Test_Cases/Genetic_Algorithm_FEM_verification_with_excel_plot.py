import numpy as np
import matplotlib.pyplot as plt


stress=np.array([-23.68421053,-57.22194764,-10.05524904,111.60285904,-119.21958634,0.,41.5674939,-4.7922403,6.26588587,-63.86575244])
strains=np.array([-7.89473684e-07,-1.90739825e-06,-3.35174968e-07,3.72009530e-06,-3.97398621e-06,0,1.38558313e-06,-1.5974134e-07,2.08862862e-07,-2.12885841e-06])

stress_from_excel=np.array([-2.37E+01,-5.72E+01,-1.01E+01,1.12E+02,-1.19E+02,0,4.16E+01,-4.79E+00,6.27E+00,-6.39E+01])
strain_from_excel=np.array([-7.89474E-07,-1.9074E-06,-3.35175E-07,3.7201E-06,-3.97399E-06,0,1.38558E-06,-1.59741E-07,2.08863E-07,-2.12886E-06])

plt.plot(stress,strains,c='b',label='From FEM')
plt.grid()
plt.plot(stress_from_excel,strain_from_excel,c='r',label='From Excel')
plt.xlabel('Stress in element')
plt.ylabel('Strains in element')
plt.title('Stress Strain plot')
plt.savefig('FEM verfication with excel and program.jpg')