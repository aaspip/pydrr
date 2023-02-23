#  DEMO script (python version) for localized DRR (DRR in windows)
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

## please first download the data from https://github.com/aaspip/data/blob/main/drr_real3d.mat
## 
#This real dataset was used previously in 
#[1] Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
#[2] Wang et al., 2021, Nonstationary predictive filtering for seismic random noise suppression — A tutorial, Geophysics, 86, W21–W30.


#load data
import scipy
from scipy import io
datas=scipy.io.loadmat("drr_real3d.mat")
dn=datas['data'];
dn=dn/dn.max();

print(np.std(dn))
print(dn.shape)

[n1,n2,n3]=dn.shape;

d1=pd.drr3d(dn,0,120,0.004,20,4);	#RR
noi1=dn-d1;

n1win=128;n2win=64;n3win=16;r1=0.5;r2=0.5;r3=0.5;
d2=pd.drr3d_win(dn,0,120,0.004,5,4,0,n1win,n2win,n3win,r1,r2,r3); #Windowed DRR or LDRR
noi2=dn-d2;

## plotting
fig = plt.figure(figsize=(10, 7))
ax=fig.add_subplot(3, 2, 1)
plt.imshow(dn.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noisy data');
ax=fig.add_subplot(3, 2, 3)
plt.imshow(d1.reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Denoised (DRR)');
ax=fig.add_subplot(3, 2, 4)
plt.imshow(noi1.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noise (DRR)');
ax=fig.add_subplot(3, 2, 5)
plt.imshow(d2.reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Denoised (LDRR)');
ax=fig.add_subplot(3, 2, 6)
plt.imshow(noi2.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noise (LDRR)');
plt.savefig('test_pydrr_drr3d_win.png',format='png',dpi=300);
plt.show()




