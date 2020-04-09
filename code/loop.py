# Secuencia de ejecucion
# f95 -O3 read_5.f metap.f iter_2.f dranxor.f 
# f95 -O3 read_fit_log.f -o b.out
# ipython loop.py
#
# El output es
# factor_propocionalidad_r Xi2_r factor_propocionalidad_i Xi2_i

import os
import numpy as np

p1min=2
p1max=3
p1dif=1

p2min=4
p2max=5
p2dif=1

p3min=15
p3max=16
p3dif=1

p4min=.3
p4max=.301
p4dif=0.1

p5min=0.1
p5max=.1001
p5dif=.05

p6min=0.
p6max=.1001
p6dif=.01

p7min=0.0
p7max=.01
p7dif=.1

v1=np.arange(p1min,p1max,p1dif)
v2=np.arange(p2min,p2max,p2dif)
v3=np.arange(p3min,p3max,p3dif)
v4=np.arange(p4min,p4max,p4dif)
v5=np.arange(p5min,p5max,p5dif)
v6=np.arange(p6min,p6max,p6dif)
v7=np.arange(p7min,p7max,p7dif)

for param1 in v1:
 for param2 in v2:
  for param3 in v3:
   for param4 in v4:
    for param5 in v5:
     for param6 in v6:
      for param7 in v7:
	print(param1,param2,param3,param4,param5,param6,param7)
	os.system('./a.out '+str(param1)+" "+str(param2)+" "+str(param3)+" "+str(param4)+" "+str(param5)+" "+str(param6)+" "+str(param7))
	os.system('./b.out ')
	os.system('mv fort.9 Output/time_av_'+str(param1)+"_"+str(param2)+"_"+str(param3)+"_"+str(param4)+"_"+str(param5)+"_"+str(param6)+"_"+str(param7)+".dat")
	os.system('mv fort.12 Output/time_all_'+str(param1)+"_"+str(param2)+"_"+str(param3)+"_"+str(param4)+"_"+str(param5)+"_"+str(param6)+"_"+str(param7)+".dat")
	os.system('mv fort.13 Output/inf_time_all_'+str(param1)+"_"+str(param2)+"_"+str(param3)+"_"+str(param4)+"_"+str(param5)+"_"+str(param6)+"_"+str(param7)+".dat")
