import os
import numpy as np

p1min=0
p1max=0.2
p1dif=0.1
p2min=0
p2max=0.6
p2dif=0.5
v1=np.arange(p1min,p1max,p1dif)
v2=np.arange(p2min,p2max,p2dif)
for param1 in v1:
 for param2 in v2:
   print(param1,param2)
   os.system('ipython reduce.py '+str(param1)+" "+str(param2)) ###Include the call to the executable, passing param 1 and param 2
