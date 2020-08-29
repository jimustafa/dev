#Modtran Plots
import matplotlib.pyplot as plt
import numpy as np
import math


file1 = open('ModtranDataTransposed.csv', 'r') 
MTD = np.loadtxt('ModtranDataTransposed.csv', delimiter=',')
file1.close()
#MTD-T= np.transpose(MTD)  

  
# plot 
plt.figure(figsize=(10,6)) # 10 is width, 6 is height
plt.plot(MTD[1],'r')
plt.plot(MTD[2],'b')
plt.plot(MTD[3],'g')
plt.title('Modtran Results')  
plt.xlabel('Wavelenght')
plt.ylabel('Transmission')
#plt.xlim(0, 15)
plt.ylim(0, 1)
plt.grid(True)
#plt.yscale('log')

plt.show()   

