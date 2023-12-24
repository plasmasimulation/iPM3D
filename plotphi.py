import numpy as np
from matplotlib import pyplot as plt
a=[]
c=[]
for i in range(0,8): 
   
    with open('./run/solve/petscphi'+str(i)+'.txt')as f:
        b=f.readline().strip()
        c.append([int(x) for x in b.split(' ')])
    

coord=np.array(c)
x=coord[0,0]+coord[1,0]
y=coord[0,1]+coord[2,1]
z=coord[0,2]+coord[4,2]
f=np.ones((x,y,z))

for i in range(0,8):
    a=np.loadtxt('./run/solve/petscphi'+str(i)+'.txt',skiprows=1)
    a=a.reshape(coord[i,0],coord[i,1],coord[i,2])
   #  print(a)

    if (i==0):
        for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x,y,z]=a[x,y,z] 
    if (i==1):
       for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x+coord[0,0],y,z]=a[x,y,z] 
    if (i==2):
      for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x,y+coord[0,1],z]=a[x,y,z]  
    if (i==3):
          for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x+coord[0,0],y+coord[0,1],z]=a[x,y,z]
    if (i==4):
      for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x,y,z+coord[0,2]]=a[x,y,z]
    if (i==5):
      for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x+coord[0,0],y,z+coord[0,2]]=a[x,y,z] 
                
    if (i==6):
      for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x,y+coord[0,1],z+coord[0,2]]=a[x,y,z]
    if (i==7):
       for x in range(0,coord[i,0]):
            for y in range(0,coord[i,1]):
               for z in range(0,coord[i,2]):
                  f[x+coord[0,0],y+coord[0,1],z+coord[0,2]]=a[x,y,z]
# print(c[0,0])
#绘制沿z方向的切片
print(f.shape)
# print(f)

xx=np.arange(0,5)
yy=np.arange(0,5)
X1,Y1=np.meshgrid(xx,yy)
# fig = plt.figure()  #定义新的三维坐标轴
fig,ax = plt.subplots(figsize=(6,8))
h=f[:,57,:]
# print(h)
sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=0, vmax=1))  
sm.set_array([])

plt.imshow(h, cmap='RdBu', norm=plt.Normalize(vmin=0, vmax=1000),extent=[0, 95, 0, 95])  
plt.colorbar(label="phi")  

plt.xlabel('X')  
plt.ylabel('Y')  
plt.title('phi x-y surface z=23')  
print(h)
plt.savefig("kk.png")

k =np.zeros((25,5))

# np.savetxt("/home/sapphire/sparta/iPM3D/testpetsc/input/geometry.txt",k,fmt='%f',delimiter=',')
                 
  
# print(b)
#打印出x方向的电场，除去边界值外，内部为一匀强电场