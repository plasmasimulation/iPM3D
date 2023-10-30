import numpy as np
a=[]
b=np.ones((6,6,6))
for i in range(0,8): 
    a=np.loadtxt('./run/ext/field/fieldx'+str(i)+'.txt')
   
    a=a.reshape(4,4,4)

    if (i==0):
        for x in range(0,3):
            for y in range(0,3):
               for z in range(0,3):
                  b[x,y,z]=a[z,y,x] 
    if (i==1):
        for x in range(3,6):
            for y in range(0,3):
               for z in range(0,3):
                  b[x,y,z]=a[z,y,x-2]
    if (i==2):
        for x in range(0,3):
            for y in range(3,6):
               for z in range(0,3):
                  b[x,y,z]=a[z,y-2,x]  
    if (i==3):
        for x in range(3,6):
            for y in range(3,6):
               for z in range(0,3):
                  b[x,y,z]=a[z,y-2,x-2] 
    if (i==4):
        for x in range(0,3):
            for y in range(0,3):
               for z in range(3,6):
                  b[x,y,z]=a[z-2,y,x]
       
    if (i==5):
        for x in range(3,6):
            for y in range(0,3):
               for z in range(3,6):
                  b[x,y,z]=a[z-2,y,x-2]
                
    if (i==6):
        for x in range(0,3):
            for y in range(3,6):
               for z in range(3,6):
                  b[x,y,z]=a[z-2,y-2,x]  
    if (i==7):
        for x in range(3,6):
            for y in range(3,6):
               for z in range(3,6):
                  b[x,y,z]=a[z-2,y-2,x-2] 
                 
  
print(b)
#打印出x方向的电场，除去边界值外，内部为一匀强电场