import numpy as np
a=[]
c=[]
for i in range(0,8): 
   
    with open('./run/logfile'+str(i)+'.txt')as f:
        b=f.readline().strip()
        c.append([int(x) for x in b.split(' ')])
    

coord=np.array(c)
x=coord[0,0]+coord[1,0]
y=coord[0,1]+coord[2,1]
z=coord[0,2]+coord[4,2]
f=np.ones((x,y,z))

for i in range(0,8):
    a=np.loadtxt('./run/logfile'+str(i)+'.txt',skiprows=1)
    a=a.reshape(coord[i,0],coord[i,1],coord[i,2])

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
                  f[x,y,z]=a[x,y,z]
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
print(f)



                 
  
# print(b)
#打印出x方向的电场，除去边界值外，内部为一匀强电场