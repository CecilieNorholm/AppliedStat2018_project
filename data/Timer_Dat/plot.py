import numpy as np
from matplotlib.pyplot import imshow, show 
import pylab as plt 

data = np.genfromtxt("/Users/Runekjaersgaard/Desktop/AppStat/Projekt/Timer_Dat/timer_Z1.dat")

l,flux = data[:,0], data[:,1]


#plt.axis([3620, 3630, -2.5e-16, 1e-15])

#plt.show()

plt.xlabel('Timestep')
plt.ylabel('Time')
plt.title('Pendulum timing')
plt.scatter(l,flux)
plt.show()
