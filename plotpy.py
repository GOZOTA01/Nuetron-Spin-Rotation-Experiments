import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px

mydata = pd.read_csv("testdata.dat", header=0)
mydata3 = np.genfromtxt('testdata.dat',delimiter=',')
data = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])

for i in range(0,len(mydata3)-1):
	x = int(mydata3[i][0])
	y = int(mydata3[i][1])
	data[x][y] = mydata3[i][2] 
	
mydata2 = mydata.to_numpy(dtype=float)
#print(mydata3)
fig = px.imshow(data)
fig.show()

