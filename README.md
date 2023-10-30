# iPM3D
plasma simulation 3d   
schedule   
test domain 3d

test domain 使用方法：  
cd test  
sh  ./setup.sh  
之后   
cd build    
mpirun -np 8 ./executablefile

test particlecmm 使用方法：  
cd testparticle  
sh  ./setup.sh  
之后  
cd build  
mpirun -np 8 ./particlecomm  

test field使用方法: 
cd testfield  
sh ./setup.sh
之后
cd run 
mpirun -np 8 ./fieldcomm
  
然后回到testfield目录   
python3 fieldplot.py  
会打印出整个区域x方向的电场

