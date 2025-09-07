import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import subprocess
dir_theta = subprocess.run(["pwd"],capture_output=True,text=True)
str_dir_theta = str(dir_theta.stdout).strip()
print(str_dir_theta+"p")

print("filename1 = ?")
filename1 = str(str_dir_theta+"/"+str(input()).strip())
print("p"+filename1)
filename1 = "".join(filename1)
print(filename1)

print("filename2 = ?")
filename2= str(str_dir_theta+"/"+str(input()).strip())
filename2="".join(filename2)
print(filename2)

print("filename3 = ?")
filename3 = str(str_dir_theta+"/"+str(input()).strip())
filename3 ="".join(filename3)
print(filename3)

#mydata = pd.read_csv("testdata2.dat", header=0)
#myxy = pd.read_csv("testdataxy.dat", header=0)
mydata = np.loadtxt(filename1, skiprows=1)
myxy = np.loadtxt("xygrid.txt", skiprows=0)

#print(mydata)
#mydata2 = mydata.to_numpy(dtype=float)

#print(myxy[:,0])
#print(myxy[:,1])
#fig = px.imshow(mydata)
fig = px.imshow(mydata,x=myxy[:,0],y=myxy[:,1],labels = dict(
x = "x(cm)",
y = "y(cm)",
color = "thetax"))
fig.show()

#filename2 = input()
mydata = np.loadtxt(filename2, skiprows=1)
myxy = np.loadtxt("xygrid.txt", skiprows=0)

#print(mydata)
#mydata2 = mydata.to_numpy(dtype=float)

#print(myxy[:,0])
#print(myxy[:,1])
#fig = px.imshow(mydata)
fig = px.imshow(mydata,x=myxy[:,0],y=myxy[:,1],labels = dict(
x = "x(cm)",
y = "y(cm)",
color = "thetay"))
fig.show()

#filename3 = input()
mydata = np.loadtxt(filename3, skiprows=1)
myxy = np.loadtxt("xygrid.txt", skiprows=0)

#print(mydata)
#print('myxy ',myxy[:,0])
#print('myxy ',myxy[:,1])

#mydata2 = mydata.to_numpy(dtype=float)

#print(myxy[:,0])
#print(myxy[:,1])
#fig = px.imshow(mydata)
fig = px.imshow(mydata,x=myxy[:,0],y=myxy[:,1],labels = dict(
x = "x(cm)",
y = "y(cm)",
color = "thetaz"))
fig.show()
