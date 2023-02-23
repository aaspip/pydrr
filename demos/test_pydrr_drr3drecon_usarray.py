#  Demo script for teleseismic denoising and reconstruction
#  as introduced in Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
#  Another version of this example is at 
#  https://github.com/chenyk1990/reproducible_research/tree/master/nonMada/usarray
# 
#  This takes about 10 minutes
#
#  Written by Yangkang Chen
#  Feb, 2018
#  Modified on Dec, 2020
#  Further polished on July, 2022
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
# REFERENCES
#  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
#  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
#  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
#  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
#  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
#  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
#  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
#  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.

import numpy as np
import matplotlib.pyplot as plt
import pydrr as pd #pd: DRR

## Please download data from https://github.com/aaspip/data/blob/main/usarray_200901181411_wfm.mat
## 
#This real dataset was used previously in 
#  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.


#load data
import scipy
from scipy import io
datas=scipy.io.loadmat("usarray_200901181411_wfm.mat")
# d, dists(shot receiver distance/offset in degree), stla, stlo, t 
d=datas['d'];
stlo=datas['stlo'];
stla=datas['stla']

## rm bad trace
inds=[17,40,69];
d=np.delete(d,inds,1);
d=pd.scale(d);
stlo=np.delete(stlo,inds,0);
stla=np.delete(stla,inds,0);
d0=d[:,104:433];
stlo0=stlo[104:433];
stla0=stla[104:433];

## rm bad trace
# inds=[18,41,70];
# d(:,inds)=[];d=drr_scale(d);
# 
# stlo(inds)=[];
# stla(inds)=[];
# d0=d(:,105:433);
# stlo0=stlo(105:433);
# stla0=stla(105:433);

## 3D processing/reconstruction
mla=[33,49];
mlo=[-116,-102];
#binning
d3d,x1,y1,mask=pd.bin3d(d0,stlo0,stla0,16,28,mlo[0],mla[0],mlo[1],mla[1]);
[stlo1,stla1]=np.meshgrid(x1,y1);


# figure;plot(stlo0,stla0,'bv');hold on;
# plot(stlo1(:),stla1(:),'r*');

# figure;imagesc(squeeze(mask(1,:,:))');colorbar;set(gca,'YDir','normal');

# figure;
# subplot(1,2,1);imagesc(squeeze(d3d(:,13,:)));caxis([-0.05,0.05]);colormap(gray);
# subplot(1,2,2);imagesc(squeeze(d3d(:,:,10)));caxis([-0.05,0.05]);colormap(gray);


## test if mask is correct
tt=(d3d*mask-d3d);
np.linalg.norm(tt.flatten()) #0->correct

ratio=np.where(mask.flatten()==1)[0].size/mask.size;
print('Sampling ratio is ',ratio);

# figure;
# subplot(1,2,1);imagesc(squeeze(mask(:,13,:)));caxis([0,1]);colormap(jet);colorbar;
# subplot(1,2,2);imagesc(squeeze(mask(:,:,10)));caxis([0,1]);colormap(jet);colorbar;

## global processing
[n1,n2,n3]=d3d.shape;
flow=0;fhigh=0.5;dt=1;N=8;Niter=10;mode=1;verb=1;eps=0.00001;K=4;
a=np.linspace(1,0,10); #linearly decreasing
d1=pd.drr3drecon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,mode,a,verb);

## local processing (first way, directly utilize drr_win3dmask.m; second way, use drr3drecon_win.m, which is an integrated function, instead)
# [n1,n2,n3]=d3d.shape;
# param.dt=dt;
# param.flow=flow;
# param.fhigh=fhigh;
# param.N=N;
# param.K=4;
# param.niter=Niter;
# param.a=a; #linearly decreasing;
# param.eps=0.00001;
# param.mode=1;
# param.verb=1;
# n1win=2000;n2win=n2;n3win=n3;
# r1=0.5;r2=0.5;r3=0.5;
# d2=drr_win3dmask(@localdrr3drecon, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);
# 
# param.amode=2;
# d3=pd.win3dmask(@localdrr3drecon_auto, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);

## fixed-rank
[n1,n2,n3]=d3d.shape;
n1win=2000;n2win=n2;n3win=n3;r1=0.5;r2=0.5;r3=0.5;
d2=pd.drr3drecon_win(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,mode,a,verb,n1win,n2win,n3win,r1,r2,r3);

## automatically chosen ranks (cleaner)
amode=2;eps2=0;#if amode=2, then eps3 is useless
#N=[4,N]; #N can also have a lower limit
d3=pd.drr3drecon_win_auto(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,mode,a,verb,amode,eps2,n1win,n2win,n3win,r1,r2,r3);

## plotting
ilon=12;
fig = plt.figure(figsize=(8, 7))
ax=fig.add_subplot(4, 1, 1)
plt.imshow(np.squeeze(d3d[:,ilon,:]),cmap='gray',clim=(-0.02, 0.02),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Raw data');
ax=fig.add_subplot(4, 1, 2)
plt.imshow(np.squeeze(d1[:,ilon,:]),cmap='gray',clim=(-0.02, 0.02),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('GDRR');
ax=fig.add_subplot(4, 1, 3)
plt.imshow(np.squeeze(d2[:,ilon,:]),cmap='gray',clim=(-0.02, 0.02),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('LDRR');
ax=fig.add_subplot(4, 1, 4)
plt.imshow(np.squeeze(d3[:,ilon,:]),cmap='gray',clim=(-0.02, 0.02),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('LDRRA');

plt.savefig('test_pydrr_drr3drecon_usarray.png',format='png',dpi=300);
plt.show()

