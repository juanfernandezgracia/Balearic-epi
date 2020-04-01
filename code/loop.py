import os
import numpy as np

p1min=0
p1max=10
p1dif=1
p2min=0
p2max=20
p2dif=1
p3min=0
p3max=1
p3dif=0.1
v1=np.arange(p1min,p1max,p1dif)
v2=np.arange(p2min,p2max,p2dif)
v3=np.arange(p3min,p3max,p3dif)

for param1 in v1:
 for param2 in v2:
  for param3 in v3:
   print(param1,param2,param3)
   os.system('./a.out '+str(param1)+" "+str(param2)+" "+str(param3))
   os.system('m v fort.9 time_av_'+str(param1)+"_"+str(param2)+"_"+str(param3)+".dat")
   os.system('m v fort.12 time_all_'+str(param1)+"_"+str(param2)+"_"+str(param3)+".dat")
