#  DEMO script (python version) for DRR
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

## generate the synthetic data
a1=np.zeros([300,80])
[n,m]=a1.shape
a3=a1.copy();
a4=a1.copy();

k=-1;
a=0.1;
b=1;
pi=np.pi

ts=np.arange(-0.055,0.055+0.002,0.002)
b1=np.zeros([len(ts)])
b2=np.zeros([len(ts)])
b3=np.zeros([len(ts)])
b4=np.zeros([len(ts)])

for t in ts:
    k=k+1;
    b1[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));
    b2[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b3[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b4[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));

t1=np.zeros([m],dtype='int')
t3=np.zeros([m],dtype='int')
t4=np.zeros([m],dtype='int')
for i in range(m):
  t1[i]=np.round(140);
  t3[i]=np.round(-2*i+220);
  t4[i]=np.round(2*i+10);
  a1[t1[i]:t1[i]+k+1,i]=b1; 
  a3[t3[i]:t3[i]+k+1,i]=b1; 
  a4[t4[i]:t4[i]+k+1,i]=b1; 

dc=a1[0:300,:]+a3[0:300,:]+a4[0:300,:];

## add noise
[n1,n2]=dc.shape
np.random.seed(201415)
n=0.1*np.random.randn(n1,n2);
dn=dc+n;
print(np.std(dn))

## decimate traces
ratio=0.7;
mask=pd.genmask(dn,ratio,'c',201415);
d0=dn*mask;
print(np.std(d0))

## Recon
flow=0;fhigh=125;dt=0.004;N=3;NN=3;Niter=10;mode=1;a=np.linspace(1,0,10);verb=1;eps=0.00001;
d1=pd.drr3drecon(d0,mask,flow,fhigh,dt,N,100,Niter,eps,mode,a,verb);
d2=pd.drr3drecon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,mode,a,verb);
noi1=dc-d1;
noi2=dc-d2;


print('SNR of RR is %g'%pd.snr(dc,d1));
print('SNR of DRR is %g'%pd.snr(dc,d2));

## compare with matlab
import scipy
from scipy import io
datas = {"d0":d0,"dc":dc,"mask":mask,"dn": dn, "d1": d1, "noi1": noi1, "d2":d2, "noi2":noi2}
scipy.io.savemat("datas2d.mat", datas)

## plot results
fig = plt.figure(figsize=(8, 8.5))
ax = fig.add_subplot(3,2,1)
plt.imshow(dc,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Clean data');
ax = fig.add_subplot(3,2,2)
plt.imshow(d0,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Incomplete noisy data');
ax = fig.add_subplot(3,2,3)
plt.imshow(d1,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Reconstructed (RR, SNR=%.4g dB)'%pd.snr(dc,d1));
ax = fig.add_subplot(3,2,4)
plt.imshow(noi1,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Error (RR)');
ax = fig.add_subplot(3,2,5)
plt.imshow(d2,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Reconstructed (DRR, SNR=%.4g dB)'%pd.snr(dc,d2));
ax = fig.add_subplot(3,2,6)
plt.imshow(noi2,cmap='jet',clim=(-0.1, 0.1),aspect=0.2);ax.set_xticks([]);ax.set_yticks([]);
plt.title('Error (DRR)');
plt.savefig('test_pydrr_drr2drecon.png',format='png',dpi=300);
plt.show()











