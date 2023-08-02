import numpy as np
import os
alf_info={}
alf_info['name']='14benz'
alf_info['nsubs']=[5,6]
alf_info['nblocks']=np.sum(alf_info['nsubs'])
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']=os.environ['CHARMMEXEC']
alf_info['temp']=298.15
