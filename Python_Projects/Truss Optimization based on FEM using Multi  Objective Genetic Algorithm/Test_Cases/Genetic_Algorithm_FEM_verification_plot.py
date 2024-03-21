import numpy as np
import matplotlib.pyplot as plt

displacement_at_each_node=np.arange(12)
stress_at_each_element=np.arange(11)

stress=np.array([-210.90153277,122.43180057,-44.2348661,62.55754757,-173.14471282,-88.4697322,122.43180057,-44.2348661,-173.14471282,62.55754757,-210.90153277])
displacement=np.array([0.,0.,7.14285714,-9.03863712,5.24707717,-16.2964926,5.24707717,-20.08805255,10.49415433,0.,3.35129719,-9.03863712])

displacement_from_paper=np.array([0,0,7.1429,-9.0386,5.2471,-16.2965,5.2471,-20.0881,10.4942,0,3.3513,-9.0386])
stress_from_paper=np.array([-210.9015,122.4318,-44.2349,62.5575,-173.1447,-88.4697,122.4318,-44.2349,-173.1447,62.5575,-210.9015])

fig,ax=plt.subplots(2,sharex='col')

ax[0].plot(displacement_at_each_node,displacement_from_paper,c='b', label='From paper')
ax[0].plot(displacement_at_each_node,displacement,c='r', label='From FEM')
ax[0].set_xlabel('Total number of displacements in the truss')
ax[0].set_ylabel('Displacements at each node')
ax[0].set_title('Stress and Displacement plots')
leg = ax[0].legend()

ax[1].plot(stress_at_each_element,stress_from_paper,c='b', label='From paper')
ax[1].plot(stress_at_each_element,stress,c='r', label='From FEM')
ax[1].set_xlabel('Total number of elements in the truss')
ax[1].set_ylabel('Stress at each node')
leg = ax[1].legend()

plt.savefig('FEM verfication with paper and program.jpg')