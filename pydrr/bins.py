def bin3d(din,x,y,nx,ny,ox=None,oy=None,mx=None,my=None):
	"""
	bin3d (previous cykbin3d or yc_bin3d): 3D seismic data binning (including 1D row vector and 2D seismics)
	IN    d:   	intput 2D data
         x:     input x coordinates
         y:     input y coordinates
         nx:    input number of binned x points
         ny:    input number of binned y points
         ox:    min of x
         oy:    min of y
         mx:    max of x
         my:    max of y
 
	 OUT   d1:  	output data
         xout:  output x coordinates
         yout:  output x coordinates
         mask:  mask operator for interpolation
 
   Copyright (C) 2015 The University of Texas at Austin
   Copyright (C) 2015 Yangkang Chen
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details: http://www.gnu.org/licenses/
 
  REFERENCES
   Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
   Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
   Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
   Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
   Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
   Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
   Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
   Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
 
   see also: yc_bin2d
	"""
	import numpy as np
	[n1,n2]=din.shape;

	if ox is None:
		ox=np.min(x);
		oy=np.min(y);
		mx=np.max(x);
		my=np.max(y);
		dx=(mx-ox)/(nx-1);
		dy=(my-oy)/(ny-1);
	else:
		dx=(mx-ox)/(nx-1);
		dy=(my-oy)/(ny-1);

	xout=np.linspace(ox,mx,nx);
	yout=np.linspace(oy,my,ny);


	dout=np.zeros([n1,nx,ny]);
	mask=np.ones([nx,ny]);

	for iy in range(0,ny):
		for ix in range(0,nx):
			inds=np.where((x>=xout[ix]) & (x<xout[ix]+dx) & (y>=yout[iy]) & (y<yout[iy]+dy))[0];
			n=len(inds);
			if n==1:
# 				print(inds.shape,dout[:,ix,iy].shape,din[:,inds].shape,np.squeeze(din[:,inds]).shape)
				dout[:,ix,iy]=np.squeeze(din[:,inds]);
			else:
				if n==0:
					mask[ix,iy]=0;
					dout[:,ix,iy]=np.zeros([n1]);
				if n>=2:
					if x[inds[0]]==xout[ix] and y[inds[0]]==yout[iy]:
						dout[:,ix,iy] = din[:,inds[0]];
					else:
						t1=np.sqrt( np.power(x[inds[0]] - xout[ix],2) + np.power(y[inds[0]] - yout[iy],2) );
						t2=np.sqrt( np.power(x[inds[1]] - xout[ix],2) + np.power(y[inds[1]] - yout[iy],2) );
						dout[:,ix,iy] = (t1*din[:,inds[1]] + t2*din[:,inds[0]])/(t1+t2);


	mask=np.multiply(np.ones([n1,1]),mask.reshape([1,nx*ny],order='F'));
	mask=mask.reshape([n1,nx,ny],order='F');


	return dout,xout,yout,mask