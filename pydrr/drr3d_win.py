def drr3d_win(D,flow=1,fhigh=124,dt=0.004,N=1,K=4,verb=0,n1win=None,n2win=None,n3win=None,r1=0.5,r2=0.5,r3=0.5):
	"""
	DRR3D_WIN: DRR3D in windows
	
	IN   D:   	intput 3D data
	     flow:   processing frequency range (lower)
	     fhigh:  processing frequency range (higher)
	     dt:     temporal sampling interval
	     N:      number of singular value to be preserved
	     K:      damping factor
	     verb:   verbosity flag (default: 0)
	     n1win:  window size
	     n2win:  window size
	     n3win:  window size
	     r1:  	overlapping ratio (default, 0.5)
	     r2:  	overlapping ratio (default, 0.5)
	     r3:  	overlapping ratio (default, 0.5)
	
	OUT  D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	REFERENCES
	Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
	Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
	Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
	Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	 
	DEMO: demos/test_pydrr_drr3d_win.m
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)
	print('n1win=',n1win,'n2win=',n2win,'r1=',r1,'r2=',r2,'r3=',r3)
	
	from .window import win3d
	from .localfun import localdrr3d
	import numpy as np
	
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	
	[nt,nx,ny]=D.shape
	print(D.shape)
	if n1win is None or n2win is None or n3win is None:
		n1win=nt;
		n2win=nx;
		n3win=ny;
		r1=0.5;
		r2=0.5;
		r3=0.5;
	param={'dt':dt,'flow':flow,'fhigh':fhigh,'N':N,'K':K,'verb':verb}
	#The N and K follows the definitions in Bai et al. (2020)
	
	D1=win3d(localdrr3d, param, D, n1win, n2win, n3win, r1, r2, r3);
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	return D1

def drr3d_win_auto(D,flow=1,fhigh=124,dt=0.004,N=1,K=4,verb=0,mode=2,n1win=None,n2win=None,n3win=None,r1=0.5,r2=0.5,r3=0.5):
	"""
	DRR3D_WIN_AUTO: DRR3D in windows with automatically selected ranks
	
	IN   D:   	intput 3D data
	     flow:   processing frequency range (lower)
	     fhigh:  processing frequency range (higher)
	     dt:     temporal sampling interval
	     N:      number of singular value to be preserved
	     K:      damping factor
	     verb:   verbosity flag (default: 0)
	     mode:   amplitude difference [1] or ratio [2]
	     n1win:  window size
	     n2win:  window size
	     n3win:  window size
	     r1:  	overlapping ratio (default, 0.5)
	     r2:  	overlapping ratio (default, 0.5)
	     r3:  	overlapping ratio (default, 0.5)
	
	OUT  D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	REFERENCES
	Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
	Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
	Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
	Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	 
	DEMO: demos/test_pydrr_drr3d_win.m
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)
	print('n1win=',n1win,'n2win=',n2win,'r1=',r1,'r2=',r2,'r3=',r3)
	
	from .window import win3d
	from .localfun import localdrr3d_auto
	import numpy as np
	
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	
	[nt,nx,ny]=D.shape
	print(D.shape)
	if n1win is None or n2win is None or n3win is None:
		n1win=nt;
		n2win=nx;
		n3win=ny;
		r1=0.5;
		r2=0.5;
		r3=0.5;
	param={'dt':dt,'flow':flow,'fhigh':fhigh,'N':N,'K':K,'verb':verb,'mode':mode}
	#The N and K follows the definitions in Bai et al. (2020)
	
	D1=win3d(localdrr3d_auto, param, D, n1win, n2win, n3win, r1, r2, r3);
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	return D1
	
def drr3drecon_win(D,mask,flow=1,fhigh=124,dt=0.004,N=1,K=4,Niter=10,eps=0.00001,mode=1,a=1,verb=0,n1win=None,n2win=None,n3win=None,r1=0.5,r2=0.5,r3=0.5):
	"""
	DRR3DRECON_WIN: DRR3DRECON in windows
	
	IN   D:   	intput 3D data
	     flow:   processing frequency range (lower)
	     fhigh:  processing frequency range (higher)
	     dt:     temporal sampling interval
	     N:      number of singular value to be preserved
	     K:      damping factor
	     verb:   verbosity flag (default: 0)
	     mode:   mode=1: denoising and reconstruction
	            mode=0: reconstruction only
	     a:      weight vector
	     n1win:  window size
	     n2win:  window size
	     n3win:  window size
	     r1:  	overlapping ratio (default, 0.5)
	     r2:  	overlapping ratio (default, 0.5)
	     r3:  	overlapping ratio (default, 0.5)
	
	OUT  D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	REFERENCES
	Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
	Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
	Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
	Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	 
	DEMO: demos/test_pydrr_drr3d_win.m
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)
	print('n1win=',n1win,'n2win=',n2win,'r1=',r1,'r2=',r2,'r3=',r3)
	
	from .window import win3dmask
	from .localfun import localdrr3drecon
	
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	[nt,nx,ny]=D.shape

	if n1win is None or n2win is None or n3win is None:
		n1win=nt;
		n2win=nx;
		n3win=ny;
		r1=0.5;
		r2=0.5;
		r3=0.5;
	param={'dt':dt,'flow':flow,'fhigh':fhigh,'N':N,'K':K,'verb':verb,'mode':mode,'niter':Niter,'a':a,'eps':eps}
	#The N and K follows the definitions in Bai et al. (2020)
	
	D1=win3dmask(localdrr3drecon, mask, param, D, n1win, n2win, n3win, r1, r2, r3);
	
	return D1

def drr3drecon_win_auto(D,mask,flow=1,fhigh=124,dt=0.004,N=1,K=4,Niter=10,eps=0.00001,mode=1,a=1,verb=0,amode=2,eps2=0.001,n1win=None,n2win=None,n3win=None,r1=0.5,r2=0.5,r3=0.5):
	"""
	DRR3DRECON_WIN: DRR3DRECON in windows
	
	IN   D:   	intput 3D data
	     flow:   processing frequency range (lower)
	     fhigh:  processing frequency range (higher)
	     dt:     temporal sampling interval
	     N:      number of singular value to be preserved
	     K:      damping factor
	     verb:   verbosity flag (default: 0)
	     mode:   mode=1: denoising and reconstruction
	            mode=0: reconstruction only
	     a:      weight vector
	     n1win:  window size
	     n2win:  window size
	     n3win:  window size
	     r1:  	overlapping ratio (default, 0.5)
	     r2:  	overlapping ratio (default, 0.5)
	     r3:  	overlapping ratio (default, 0.5)
	
	OUT  D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	REFERENCES
	Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
	Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
	Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
	Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	 
	DEMO: demos/test_pydrr_drr3d_win.m
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)
	print('n1win=',n1win,'n2win=',n2win,'r1=',r1,'r2=',r2,'r3=',r3)
	
	from .window import win3dmask
	from .localfun import localdrr3drecon_auto
	
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	[nt,nx,ny]=D.shape

	if n1win is None or n2win is None or n3win is None:
		n1win=nt;
		n2win=nx;
		n3win=ny;
		r1=0.5;
		r2=0.5;
		r3=0.5;
	param={'dt':dt,'flow':flow,'fhigh':fhigh,'N':N,'K':K,'verb':verb,'mode':mode,'amode':amode,'niter':Niter,'a':a,'eps':eps,'eps2':eps2}
	#The N and K follows the definitions in Bai et al. (2020)
	
	D1=win3dmask(localdrr3drecon_auto, mask, param, D, n1win, n2win, n3win, r1, r2, r3);
	
	return D1

