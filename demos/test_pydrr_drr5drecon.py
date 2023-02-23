#  DEMO script (python version) for DRR 5D reconstruction (it is also verified to work for 4D, we just not provide a 4D example)
#  
#  Copyright (C) 2022 Yangkang Chen
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details: http://www.gnu.org/licenses/
#  
## generate synthetic data
#This synthetic data was used in Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
import numpy as np
import matplotlib.pyplot as plt
import pydrr as pd #pd: DRR

#DATAPATH:
#https://github.com/chenyk1990/reproducible_research/blob/master/odrr5d/matfun/yc_synth5d.mat

#load data
import scipy
from scipy import io
datas=scipy.io.loadmat("yc_synth5d.mat")
d=datas['data5d'];
d=d/d.max();

#add noise
[n1,n2,n3,n4,n5]=d.shape;
np.random.seed(201415)
n=0.1*np.random.randn(n1,n2,n3,n4,n5);
dn=d+n;
print(np.std(dn))

## decimate traces
ratio=0.3;
mask=pd.genmask(dn.reshape(n1,n2*n3*n4*n5,order='F'),ratio,'c',201415);
mask=mask.reshape(n1,n2,n3,n4,n5,order='F');
d0=dn*mask;
print(np.std(d0))

## Recon
flow=5;fhigh=100;dt=0.004;N=4;NN=3;Niter=10;mode=1;a=np.linspace(1,0,10);verb=1;eps=0.00000000001;
d1=pd.drr5drecon(d0,mask,flow,fhigh,dt,N,100,Niter,eps,mode,a,verb);
d2=pd.drr5drecon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,mode,a,verb);
noi1=d-d1;
noi2=d-d2;

## compare SNR
snr1=pd.snr(d[:,:,:,5,5],d1[:,:,:,5,5],2);
snr2=pd.snr(d[:,:,:,5,5],d2[:,:,:,5,5],2);
print('SNR of RR is %g'%snr1);
print('SNR of DRR is %g'%snr2);

## compare with matlab
import scipy
from scipy import io
datas = {"d0":d0,"dc":d,"mask":mask,"dn": dn, "d1": d1, "noi1": noi1, "d2":d2, "noi2":noi2}
scipy.io.savemat("datas5d.mat", datas)


## plotting
fig = plt.figure(figsize=(8, 8))
ax=fig.add_subplot(3, 2, 1)
plt.imshow(dn[:,:,:,5,5].transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noisy data');
ax=fig.add_subplot(3, 2, 2)
plt.imshow(d0[:,:,:,5,5].transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Incomplete data');
ax=fig.add_subplot(3, 2, 3)
plt.imshow(d1[:,:,:,5,5].reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Reconstructed (RR, SNR=%.4g dB)'%snr1);
ax=fig.add_subplot(3, 2, 4)
plt.imshow(noi1[:,:,:,5,5].transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Error (RR)');
ax=fig.add_subplot(3, 2, 5)
plt.imshow(d2[:,:,:,5,5].reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Reconstructed (DRR, SNR=%.4g dB)'%snr2);
ax=fig.add_subplot(3, 2, 6)
plt.imshow(noi2[:,:,:,5,5].transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.05,0.05),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Error (DRR)');
plt.savefig('test_pydrr_drr5drecon.png',format='png',dpi=300);
plt.show()







