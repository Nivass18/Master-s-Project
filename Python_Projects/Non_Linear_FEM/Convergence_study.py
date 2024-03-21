import numpy as np 
import matplotlib.pyplot as plt
time=np.array([ 0.1,0.05,0.02,0.01])
no_of_elements=np.array([36,45,82,141])

plt.plot(no_of_elements,time)
plt.title('Time_step and No of Elements')
plt.ylabel('Times')
plt.xlabel('No_of_elements')
plt.savefig('Times_and_no_of_elements.png')