import numpy as np
import scipy
from scipy import linalg
def drr3d(D, flow=1, fhigh=124, dt=0.004, N=1, K=3, verb=0):
	"""
	DRR3D: 3D damped rank-reduction method for denoising (also known as FXYDMSSA)
	
	INPUT
	D:   	intput 3D data (ndarray)
	flow:   processing frequency range (lower)
	fhigh:  processing frequency range (higher)
	dt:     temporal sampling interval
	N:      number of singular value to be preserved
	K:      damping factor
	verb:   verbosity flag (default: 0)
	
	OUTPUT  
	D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	Ported to Python in 2022 by Yangkang Chen (Verified to be correct, the same as Matlab version)
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References   
	[1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	[2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	[3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	[4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	[5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	
	 DEMO
	 demos/test_pydrr_drr2d.py
	 demos/test_pydrr_drr3d.py
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)

	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)

	[nt,nx,ny]=D.shape
	D1=np.zeros([nt,nx,ny])
	
	nf=nextpow2(nt);
	nf=int(nf)
	
	#Transform into F-X domain
	DATA_FX=np.fft.fft(D,nf,0);
	DATA_FX0=np.zeros([nf,nx,ny],dtype=np.complex_);

	#First and last nts of the DFT.
	ilow  = np.floor(flow*dt*nf)+1;
	if ilow<1:
		ilow=1
	
	ihigh = np.floor(fhigh*dt*nf)+1;
	if ihigh > np.floor(nf/2)+1:
		ihigh=np.floor(nf/2)+1
	
	ilow=int(ilow)
	ihigh=int(ihigh)
	
	lx=int(np.floor(nx/2)+1);
	lxx=nx-lx+1;
	ly=int(np.floor(ny/2)+1);
	lyy=ny-ly+1;
	M=np.zeros([lx*ly,lxx*lyy],dtype=np.complex_);
	
	#main loop
	for k in range(ilow,ihigh+1):
		
		M=P_H(DATA_FX[k-1,:,:],lx,ly); 
		M=P_RD(M,N,K);
		DATA_FX0[k-1,:,:]=P_A(M,nx,ny,lx,ly);
		
		if np.mod(k,5)==0 and verb==1 : 
			print('F %d is done!\n\n'%k);
	
	for k in range(int(nf/2)+2,nf+1):
		DATA_FX0[k-1,:,:] = np.conj(DATA_FX0[nf-k+1,:,:]);

	#Back to TX (the output)
	D1=np.real(np.fft.ifft(DATA_FX0,nf,0));
	D1=D1[0:nt,:,:];
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	
	return D1

def drr3d_auto(D, flow=1, fhigh=124, dt=0.004, N=1, K=3, verb=0, mode=2, eps=0.01):
	"""
	DRR3D_AUTO: 3D damped rank-reduction method for denoising with automatically selected ranks
	
	INPUT
	D:   	intput 3D data (ndarray)
	flow:   processing frequency range (lower)
	fhigh:  processing frequency range (higher)
	dt:     temporal sampling interval
	N:      number of singular value to be preserved
	K:      damping factor
	verb:   verbosity flag (default: 0)
	mode:   1 (amplitude difference)
			2 (ratio)
	eps: threshold (fixed in default)
	
	OUTPUT  
	D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	Ported to Python in 2022 by Yangkang Chen (Verified to be correct, the same as Matlab version)
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References   
	[1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	[2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	[3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	[4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	[5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	
	 DEMO
	 demos/test_pydrr_drr2d_diffraction.py
	 demos/test_pydrr_drr3d_diffraction.py
	"""
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'verb=',verb)

	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)

	if type(N) != list:
		Nmax=N;
		Nmin=1;
	else:
		Nmin=N[0];
		Nmax=N[1];	


	[nt,nx,ny]=D.shape
	D1=np.zeros([nt,nx,ny])
	
	nf=nextpow2(nt);
	nf=int(nf)
	
	#Transform into F-X domain
	DATA_FX=np.fft.fft(D,nf,0);
	DATA_FX0=np.zeros([nf,nx,ny],dtype=np.complex_);

	#First and last nts of the DFT.
	ilow  = np.floor(flow*dt*nf)+1;
	if ilow<1:
		ilow=1
	
	ihigh = np.floor(fhigh*dt*nf)+1;
	if ihigh > np.floor(nf/2)+1:
		ihigh=np.floor(nf/2)+1
	
	ilow=int(ilow)
	ihigh=int(ihigh)
	
	lx=int(np.floor(nx/2)+1);
	lxx=nx-lx+1;
	ly=int(np.floor(ny/2)+1);
	lyy=ny-ly+1;
	M=np.zeros([lx*ly,lxx*lyy],dtype=np.complex_);
	
	#main loop
	for k in range(ilow,ihigh+1):
		
		M=P_H(DATA_FX[k-1,:,:],lx,ly); 
		M=P_RDA(M,Nmin,Nmax,K,mode,eps);
		DATA_FX0[k-1,:,:]=P_A(M,nx,ny,lx,ly);
		
		if np.mod(k,5)==0 and verb==1 : 
			print('F %d is done!\n\n'%k);
	
	for k in range(int(nf/2)+2,nf+1):
		DATA_FX0[k-1,:,:] = np.conj(DATA_FX0[nf-k+1,:,:]);

	#Back to TX (the output)
	D1=np.real(np.fft.ifft(DATA_FX0,nf,0));
	D1=D1[0:nt,:,:];
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	
	return D1


def drr3drecon(D, MASK, flow=1, fhigh=124, dt=0.004, N=3, K=3, Niter=10,eps=0.00001,mode=0,a=1,verb=0):
	"""
	DRR3DRECON: 3D damped rank-reduction method for denoising and reconstruction 
	
	INPUT  
	D:   	 intput 3D data (ndarray)
	MASK:	 sampling mask
	flow:   processing frequency range (lower)
	fhigh:  processing frequency range (higher)
	dt:     temporal sampling interval
	N:      number of singular value to be preserved
	K:      damping factor
	Niter:  number of iterations
	eps:    termination parameter
	Niter:  number of maximum iteration
	eps:    tolerence (||S(n)-S(n-1)||_F<eps*S(n))
	mode:   mode=1: denoising and reconstruction
	        mode=0: reconstruction only
	a:      linear-decreasing scale vector
	verb:   verbosity flag (default: 0)
	
	OUTPUT  
	D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	Ported to Python in 2022 by Yangkang Chen (Verified to be correct, the same as Matlab version)
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References  
	[1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	[2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	[3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	[4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	[5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	
	 DEMO
	 demos/test_pydrr_drr2drecon.py
	 demos/test_pydrr_drr3drecon.py
	"""
	

	if mode==0:
		a=np.ones([Niter,1]);
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'Niter=',Niter,'eps=',eps,'mode=',mode,'verb=',verb)
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	if MASK.ndim==2:	#for 2D problems
		MASK=np.expand_dims(MASK, axis=2)
		
	mask=np.squeeze(MASK[0,:,:]);
	if mask.ndim==1:	#for 2D problems
		mask=np.expand_dims(mask, axis=1)
		
	[nt,nx,ny]=D.shape
	D1=np.zeros([nt,nx,ny])
	
	nf=nextpow2(nt);
	nf=int(nf)
	
	#Transform into F-X domain
	DATA_FX=np.fft.fft(D,nf,0);
	DATA_FX0=np.zeros([nf,nx,ny],dtype=np.complex_);

	#First and last nts of the DFT.
	ilow  = np.floor(flow*dt*nf)+1;
	if ilow<1:
		ilow=1
	
	ihigh = np.floor(fhigh*dt*nf)+1;
	if ihigh > np.floor(nf/2)+1:
		ihigh=np.floor(nf/2)+1
	
	ilow=int(ilow)
	ihigh=int(ihigh)
	
	lx=int(np.floor(nx/2)+1);
	lxx=nx-lx+1;
	ly=int(np.floor(ny/2)+1);
	lyy=ny-ly+1;
	M=np.zeros([lx*ly,lxx*lyy],dtype=np.complex_);
	
	print('mask.shape',mask.shape)
	#main loop
	for k in range(ilow,ihigh+1):
		S_obs=np.squeeze(DATA_FX[k-1,:,:]);   
		if S_obs.ndim==1:	#for 2D problems
			S_obs=np.expand_dims(S_obs, axis=1)
		Sn_1=S_obs;
		
		for iter in range(Niter):
			
			M=P_H(Sn_1,lx,ly); 
			M=P_RD(M,N,K);
			Sn=P_A(M,nx,ny,lx,ly);
			
			Sn=a[iter]*S_obs+(1-a[iter])*mask*Sn+(1-mask)*Sn;
			
			if np.linalg.norm(Sn-Sn_1,'fro')<eps:
				break;
			
			Sn_1=Sn;
		DATA_FX0[k-1,:,:]=Sn;
		
		if np.mod(k,5)==0 and verb==1 : 
			print('F %d is done!'%k);
	
	for k in range(int(nf/2)+2,nf+1):
		DATA_FX0[k-1,:,:] = np.conj(DATA_FX0[nf-k+1,:,:]);

	#Back to TX (the output)
	D1=np.real(np.fft.ifft(DATA_FX0,nf,0));
	D1=D1[0:nt,:,:];
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	
	return D1

def drr3drecon_auto(D, MASK, flow=1, fhigh=124, dt=0.004, N=3, K=3, Niter=10,eps=0.00001,mode=0,a=1,verb=0,amode=2,eps2=0.01):
	"""
	DRR3DRECON: 3D damped rank-reduction method for denoising and reconstruction with automatically selected ranks
	
	INPUT  
	D:   	 intput 3D data (ndarray)
	MASK:	 sampling mask
	flow:   processing frequency range (lower)
	fhigh:  processing frequency range (higher)
	dt:     temporal sampling interval
	N:      number of singular value to be preserved
	K:      damping factor
	Niter:  number of iterations
	eps:    termination parameter
	Niter:  number of maximum iteration
	eps:    tolerence (||S(n)-S(n-1)||_F<eps*S(n))
	mode:   mode=1: denoising and reconstruction
	        mode=0: reconstruction only
	amode:  automatic rank strategy
			amode=1: (amplitude difference)
			amode=2: (ratio)
	a:      linear-decreasing scale vector
	verb:   verbosity flag (default: 0)
	eps2:   threshold for amplitude difference

	OUTPUT  
	D1:  	output data
	
	Copyright (C) 2013 The University of Texas at Austin
	Copyright (C) 2013 Yangkang Chen
	Modified 2015 by Yangkang Chen
	Ported to Python in 2022 by Yangkang Chen (Verified to be correct, the same as Matlab version)
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References  
	[1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
	[2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
	[3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
	[4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	[5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	
	 DEMO
	 demos/test_pydrr_drr3drecon_usarray.py
	"""
	

	if mode==0:
		a=np.ones([Niter,1]);
	print('flow=',flow,'fhigh=',fhigh,'dt=',dt,'N=',N,'K=',K,'Niter=',Niter,'eps=',eps,'mode=',mode,'verb=',verb)
	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
	if MASK.ndim==2:	#for 2D problems
		MASK=np.expand_dims(MASK, axis=2)
		
	mask=np.squeeze(MASK[0,:,:]);
	if mask.ndim==1:	#for 2D problems
		mask=np.expand_dims(mask, axis=1)
		
	
	if type(N) != list:
		Nmax=N;
		Nmin=1;
	else:
		Nmin=N[0];
		Nmax=N[1];	
	
	[nt,nx,ny]=D.shape
	D1=np.zeros([nt,nx,ny])
	
	nf=nextpow2(nt);
	nf=int(nf)
	
	#Transform into F-X domain
	DATA_FX=np.fft.fft(D,nf,0);
	DATA_FX0=np.zeros([nf,nx,ny],dtype=np.complex_);

	#First and last nts of the DFT.
	ilow  = np.floor(flow*dt*nf)+1;
	if ilow<1:
		ilow=1
	
	ihigh = np.floor(fhigh*dt*nf)+1;
	if ihigh > np.floor(nf/2)+1:
		ihigh=np.floor(nf/2)+1
	
	ilow=int(ilow)
	ihigh=int(ihigh)
	
	lx=int(np.floor(nx/2)+1);
	lxx=nx-lx+1;
	ly=int(np.floor(ny/2)+1);
	lyy=ny-ly+1;
	M=np.zeros([lx*ly,lxx*lyy],dtype=np.complex_);
	
	print('mask.shape',mask.shape)
	#main loop
	for k in range(ilow,ihigh+1):
		S_obs=np.squeeze(DATA_FX[k-1,:,:]);   
		if S_obs.ndim==1:	#for 2D problems
			S_obs=np.expand_dims(S_obs, axis=1)
		Sn_1=S_obs;
		
		for iter in range(Niter):
			
			M=P_H(Sn_1,lx,ly); 
			M=P_RDA(M,Nmin,Nmax,K,amode,eps2);
			Sn=P_A(M,nx,ny,lx,ly);
			
			Sn=a[iter]*S_obs+(1-a[iter])*mask*Sn+(1-mask)*Sn;
			
			if np.linalg.norm(Sn-Sn_1,'fro')<eps:
				break;
			
			Sn_1=Sn;
		DATA_FX0[k-1,:,:]=Sn;
		
		if np.mod(k,5)==0 and verb==1 : 
			print('F %d is done!'%k);
	
	for k in range(int(nf/2)+2,nf+1):
		DATA_FX0[k-1,:,:] = np.conj(DATA_FX0[nf-k+1,:,:]);

	#Back to TX (the output)
	D1=np.real(np.fft.ifft(DATA_FX0,nf,0));
	D1=D1[0:nt,:,:];
	
	if ny==1:	#for 2D problems
		D1=np.squeeze(D1)
	
	return D1


def nextpow2(N):
    """ Function for finding the next power of 2 """
    n = 1
    while n < N: n *= 2
    return n
    

def P_H(din,lx,ly):
	""" forming block Hankel matrix """
	[nx,ny]=din.shape;
	lxx=nx-lx+1;
	lyy=ny-ly+1;
	dout=np.zeros([lx*ly,lxx*lyy],dtype=np.complex_)
	
	for j in range(1,ny+1):
		r=scipy.linalg.hankel(din[0:lx,j-1],din[lx-1:nx,j-1]);
		if j<ly:
			for id in range(1,j+1):
				dout[(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,(id-1)*lxx:lxx+(id-1)*lxx] = r;
		else:
			for id in range(1,ny-j+2):
				dout[(ly-1)*lx-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx]=r;	
	return dout


def P_RD(din,N,K):
	"""Rank reduction on the block Hankel matrix"""
	[U,D,V]=scipy.linalg.svd(din)
	for j in range(1,N+1):
		D[j-1]=D[j-1]*(1-np.power(D[N],K)/(np.power(D[j-1],K)+0.000000000000001))
	dout=np.mat(U[:,0:N])*np.mat(np.diag(D[0:N]))*np.mat(V[0:N,:]);

	return dout

def P_RDA(din,Nmin,Nmax,K,mode,eps):
	"""Rank reduction on the block Hankel matrix with automatically chosen rank"""
	
	[U,D,V]=scipy.linalg.svd(din)
	
	if mode==1:
		N=Nmax;
		for i in range(0,Nmax):
			if((D[i]*D[i]-D[i+1]*D[i+1])/(D[i]*D[i])) >=eps:
				N=i+1;
				continue
		if N<Nmin:
			N=Nmin;
	else:
	#used in Chen et al., 2017, Geophysics;Chen et al., 2019, NC;Wang et al., 2020, Geophysics
		N=Nmax;
		df=np.zeros([Nmax,1]);
		for i in range(0,Nmax):
			df[i]= D[i]/(D[i+1]+0.000000000000001);
		N=np.where(df==np.max(df));
		
		N=N[0][0]+1;
		if N<Nmin:
			N=Nmin;
	
	for j in range(1,N+1):
		D[j-1]=D[j-1]*(1-np.power(D[N],K)/(np.power(D[j-1],K)+0.000000000000001))
	dout=np.mat(U[:,0:N])*np.mat(np.diag(D[0:N]))*np.mat(V[0:N,:]);

	return dout
	
def P_A(din,nx,ny,lx,ly):
	""" Averaging the block Hankel matrix to output the result """ 
	lxx=nx-lx+1;
	lyy=ny-ly+1;

	dout=np.zeros([nx,ny],dtype=np.complex_);
	
	for j in range(1,ny+1):
		if j<ly:
			for id in range(1,j+1):
				dout[:,j-1] =dout[:,j-1]+ ave_antid(din[(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,(id-1)*lxx:lxx+(id-1)*lxx])/j;
		else:
			for id in range(1,ny-j+2):
				dout[:,j-1] =dout[:,j-1]+ ave_antid(din[(ly-1)*lx-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx])/(ny-j+1);
	
	return dout
    
def ave_antid(din):
	""" averaging along antidiagonals """
	[n1,n2]=din.shape;
	nout=n1+n2-1;
	dout=np.zeros(nout,dtype=np.complex_);

	for i in range(1,nout+1):
		if i<n1:
			for id in range(1,i+1):
				dout[i-1]=dout[i-1] + din[i-(id-1)-1,(id-1)]/i; 
		else:
			for id in range(1,nout+2-i):
				dout[i-1]=dout[i-1] + din[n1-(id-1)-1,(i-n1)+(id-1)]/(nout+1-i); 
	return dout




    
