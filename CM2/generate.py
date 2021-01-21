
# -*- coding: utf-8 -*-

import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from skimage import filters
from  scipy.signal import medfilt
import sys
import time
import torch
from PIL import Image, ImageDraw


current_dir=os.getcwd()
#work_case='Wasserstein'
work_dir=r'/home/mamaya/Documents/First_project/GPR_curve_ray/TImage_2/Reference_model_data/Generator'#+ '\\' + work_case
os.chdir(work_dir)


# Load TI
ti_dir=work_dir
ti_file='TI_BR_1ch.png'
ti_path= ti_dir+'/'+ti_file
ti = Image.open(ti_path)
ti=np.array(ti)
if len(ti.shape)==2: # the array is 2d, convert to a 3D array
	ti=ti.reshape((ti.shape[0],ti.shape[1],1))
ti = ti.transpose( (2,0,1) )
ti = ti / 128. - 1.0

#%% Check generated models for a given epoch
DoFiltering=False
DoThreshold=True
TriCat=True
epoch=299

gpath=work_dir+'/netG_epoch_'+str(epoch)+'.pth' 

cuda=False

from generators import G as Generator
netG = Generator(nc=1, nz=3, ngf=64, gfs=3, ngpu=1,cuda=cuda, gpath=gpath)
netG.eval()

#rn_seed=2043
#np.random.seed(rn_seed)


nz=3
zx=5
zy=3
znp = np.random.uniform(0,1, (1,nz, zx, zy))*2-1

z = torch.from_numpy(znp).float()

t0=time.time()

if cuda:
    netG.cuda()
    z = z.cuda()
    
model = netG(z)
model=model.detach().cpu().numpy()
model=0.5*(model+1)

if DoFiltering==True:
    for i in range(0,model.shape[0]):
        model[i,0,:,:,:]=medfilt(model[i,0,:,:,:], kernel_size=(3,1,1))

if DoThreshold==True:
    if TriCat:
        model[model<0.334]=0
        model[model>=0.667]=2
        model[np.where((model > 0) & (model < 2))]=1
        model=model/2.0 
    else: # Binary image
        threshold=0.5
        model[model<threshold]=0
        model[model>=threshold]=1
    
print(time.time()-t0)


m=np.array(model[:,0,:,:])

#DoSaveFig=True
#fig = plt.figure(figsize=(12,12))
#plt.subplot(2,2,1)
#plt.title('Patch of TI')
#plt.imshow(ti[0,:289,:289],)
#plt.subplot(2,2,2)
#plt.title('Realz #1')
#plt.imshow(m[1,:,:],)
#plt.subplot(2,2,3)
#plt.title('Realz #2')
#plt.imshow(m[2,:,:],cmap='gray')
#plt.subplot(2,2,4)
plt.title('Realz #3')
plt.imshow(m[0,:,:])#,cmap='gray')
plt.show()
if DoSaveFig:
    fig.savefig('BR_wass_ep299_rect5.png',dpi=300)


