def win2d(oper, param, din, n1win=None, n2win=None, r1=0.5, r2=0.5):
	"""
	Processing in 2D windows
	2D version similar to process_win
	coding strategy follows exactly usr/cm/Mfxmssa_win.c
	
	din:          input data
	oper:         operator
	param:        parameters of operator
	n1win:        first window length
	n2win:        second window length
	r1:           first overlapping ratio
	r2:           second overlapping ratio
	
	dout:     output data
	
	Author      : Yangkang Chen
	
	Date        : Jan, 2018
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as pwublished
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	see also 
	win3d(),win3dmask()
	
	Example: 
	demos/test_pydrr_drr2d_win.py
	demos/test_pydrr_drr3d_win.py
	"""
	import numpy as np
	
	[n1,n2]=din.shape;
	
	if n1win is None or n2win is None:
		n1win=n1;
		n2win=n2;
		
	nov1=(1-r1)*n1win;  nov1=int(nov1); # non-overlapping size 1
	nov2=(1-r2)*n2win;  nov2=int(nov2); # non-overlapping size 2
	ov1=r1*n1win;       ov1=int(ov1);   # overlapping size 1
	ov2=r2*n2win;       ov2=int(ov2);   # overlapping size 2

	n1pad=n1win;        # padding size 1
	nw1=1;
	while n1pad<n1:
		n1pad=n1pad+nov1;nw1=nw1+1;

	n2pad=n2win;        # padding size 2
	nw2=1;
	while n2pad<n2:
		n2pad=n2pad+nov2;nw2=nw2+1;


	D1=np.zeros([n1pad,n2pad]);D1[0:n1,0:n2]=din; # copy din into padded D1
	D2=np.zeros([n1pad,n2pad]);

	for iw2 in range(0,nw2):
		for iw1 in range(0,nw1):
			s1=iw1*nov1;s2=iw2*nov2;
			dtmp=D1[s1:s1+n1win,s2:s2+n2win];
			# uncomment this line for checking the correctness (checked 100%  correct)
			dtmp = oper(dtmp,param);
			# only valid for space-independent param
			# for reconstruction, with mask, param needs to be changed
			dtmp=win_weight2d(dtmp,iw1,iw2,nw1,nw2,n1win,n2win,ov1,ov2);
			D2[s1:s1+n1win,s2:s2+n2win]=D2[s1:s1+n1win,s2:s2+n2win]+dtmp;
	dout=D2[0:n1,0:n2];
	
	return dout


def win_weight2d(din,iw1,iw2,nw1,nw2,n1win,n2win,ov1,ov2):
	"""
	Getting the tapering weight 
	follow exactly usr/cm/win.c
	
	Author      : Yangkang Chen
	Date        :  July, 2016 (C)
				   Jan,  2018 (Matlab)
				   July, 2022 (Python)
	
	float din /*input data*/,
	int iw1 /*starting window 1 in dst*/,
	int iw2 /*starting window 2 in dst*/,
	int nw1 /*no of windows 1 in src*/,
	int nw2 /*no of windows 2 in src*/,
	int n1win /*window length 1 in src*/,
	int n2win /*window legnth 2 in src*/,
	int ov1 /*copy length in axis1*/,
	int ov2 /*copy length in axis2*/)
	"""
	
	if iw2!=0:
		for i1 in range(0,n1win):
			for i2 in range(0,ov2):
				din[i1,i2]=din[i1,i2]*(i2+1)/(ov2+1);
	if iw2!=nw2-1:
		for i1 in range(0,n1win):
			for i2 in range(0,ov2):
				din[i1,n2win-ov2+i2] = din[i1,n2win-ov2+i2]*(ov2-i2)/(ov2+1);
	if iw1!=0:
		for i2 in range(0,n2win):
			for i1 in range(0,ov1):
				din[i1,i2] = din[i1,i2]*(i1+1)/(ov1+1);
	if iw1!=nw1-1:
		for i2 in range(0,n2win):
			for i1 in range(0,ov1):
				din[n1win-ov1+i1,i2] = din[n1win-ov1+i1,i2]*(ov1-i1)/(ov1+1);
	dout=din;
	return dout

	
	
def win3d(oper, param, din, n1win=None, n2win=None, n3win=None, r1=0.5, r2=0.5, r3=0.5):
	"""
	DRR_WIN3D: Processing in 3D windows
	2D version similar to process_win
	coding strategy follows exactly usr/cm/Mfxmssa_win.c
	
	din:          input data
	oper:         operator
	param:        parameters of operator
	n1win:        first window length
	n2win:        second window length
	n3win:        third window length
	r1:           first overlapping ratio
	r2:           second overlapping ratio
	r3:           third overlapping ratio
	
	dout:     output data
	
	Author      : Yangkang Chen
	
	Date        : Jan, 2018
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	see also 
	win2d(),win3dmask()
	
	Example: 
	demos/test_pydrr_drr2d_win.py
	demos/test_pydrr_drr3d_win.py
	
	
	 REFERENCES
	 Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
	 Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
	 Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
	"""
	import numpy as np
	[n1,n2,n3]=din.shape;

	if n1win is None or n2win is None or n3win is None:
		n1win=n1;
		n2win=n2;
		n3win=n3;
		
	nov1=(1-r1)*n1win; nov1=int(nov1); # non-overlapping size 1
	nov2=(1-r2)*n2win; nov2=int(nov2); # non-overlapping size 2
	nov3=(1-r3)*n3win; nov3=int(nov3); # non-overlapping size 3
	ov1=r1*n1win;      ov1=int(ov1);   # overlapping size 1
	ov2=r2*n2win;      ov2=int(ov2);  # overlapping size 2
	ov3=r3*n3win;      ov3=int(ov3);  # overlapping size 3
	
	n1pad=n1win;        # padding size 1
	nw1=1;
	while n1pad<n1:
		n1pad=n1pad+nov1;nw1=nw1+1;

	n2pad=n2win;        # padding size 2
	nw2=1;
	while n2pad<n2:
		n2pad=n2pad+nov2;nw2=nw2+1;

	n3pad=n3win;        # padding size 3
	nw3=1;
	while n3pad<n3:
		n3pad=n3pad+nov3;nw3=nw3+1;

	D1=np.zeros([n1pad,n2pad,n3pad]);D1[0:n1,0:n2,0:n3]=din;# copy din into padded D1
	D2=np.zeros([n1pad,n2pad,n3pad]);

	for iw3 in range(0,nw3):
		for iw2 in range(0,nw2):
			for iw1 in range(0,nw1):
				s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;
				dtmp=D1[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win];
# 				print(dtmp.shape)
				# uncomment this line for checking the correctness (checked 100%  correct)
				dtmp = oper(dtmp,param);
				if dtmp.ndim==2:	#for 2D problems
					dtmp=np.expand_dims(dtmp, axis=2)
				
# 				print(dtmp.shape)
				# only valid for space-independent param
				# for reconstruction, with mask, param needs to be changed
				
				dtmp=win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
				D2[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win]=D2[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win]+dtmp;
	dout=D2[0:n1,0:n2,0:n3];

	return dout


def win_weight3d(din,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3):
	"""
	Getting the tapering weight 
	follow exactly usr/cm/win.c
	
	Author      : Yangkang Chen
	Date        :  July, 2016 (C)
				   Jan,  2018 (Matlab)
				   July, 2022 (Python)
				   
	float din /*input data*/,
	int iw1 /*starting window 1 in dst*/,
	int iw2 /*starting window 2 in dst*/,
	int iw3 /*starting window 3 in dst*/,
	int nw1 /*no of windows 1 in src*/,
	int nw2 /*no of windows 2 in src*/,
	int nw3 /*no of windows 3 in src*/,
	int n1win /*window length 1 in src*/,
	int n2win /*window legnth 2 in src*/,
	int n3win /*window legnth 3 in src*/,
	int ov1 /*copy length in axis1*/,
	int ov2 /*copy length in axis2*/,
	int ov3 /*copy length in axis3*/)

	"""
	
# 	print(din.shape)
	if iw3!=0:
		for i1 in range(0,n1win):
			for i2 in range(0,n2win):
				for i3 in range(0,ov3):
					din[i1,i2,i3] = din[i1,i2,i3]*(i3+1)/(ov3+1);

	if iw3!=nw3-1:
		for i1 in range(0,n1win):
			for i2 in range(0,n2win):
				for i3 in range(0,ov3):
					din[i1,i2,n3win-ov3+i3] = din[i1,i2,n3win-ov3+i3]*(ov3-i3)/(ov3+1);

	if iw2!=0:
		for i3 in range(0,n3win):
			for i1 in range(0,n1win):
				for i2 in range(0,ov2):
					din[i1,i2,i3] = din[i1,i2,i3]*(i2+1)/(ov2+1);

	if iw2!=nw2-1:
		for i3 in range(0,n3win):
			for i1 in range(0,n1win):
				for i2 in range(0,ov2):
					din[i1,n2win-ov2+i2,i3] = din[i1,n2win-ov2+i2,i3]*(ov2-i2)/(ov2+1);

	if iw1!=0:
		for i3 in range(0,n3win):
			for i2 in range(0,n2win):
				for i1 in range(0,ov1):
					din[i1,i2,i3] = din[i1,i2,i3]*(i1+1)/(ov1+1);

	if iw1!=nw1-1:
		for i3 in range(0,n3win):
			for i2 in range(0,n2win):
				for i1 in range(0,ov1):
					din[n1win-ov1+i1,i2,i3] = din[n1win-ov1+i1,i2,i3]*(ov1-i1)/(ov1+1);
	dout=din;
	return dout
	
	
	
def win2dmask():



	return
	
	
	
def win3dmask(oper, mask, param, din, n1win=None, n2win=None, n3win=None, r1=0.5, r2=0.5, r3=0.5):
	"""
	Processing in 3D windows (can also deal with 2D processing)
	2D version similar to process_win
	coding strategy follows exactly usr/cm/Mfxmssa_win.c
	
	din:          input data
	oper:         operator
	mask:         sampling/masking operator for data reconstruction
	param:        parameters of operator
	n1win:        first window length
	n2win:        second window length
	n3win:        third window length
	r1:           first overlapping ratio
	r2:           second overlapping ratio
	r3:           third overlapping ratio
	
	dout:     output data
	
	Author      : Yangkang Chen
	
	Date        : Jan, 2018
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	see also 
	win2d(),win3d()
	
	Example: 
	demos/test_pydrr_drr2d_win.py
	demos/test_pydrr_drr3d_win.py
	"""
	import numpy as np
	[n1,n2,n3]=din.shape;

	if n1win is None or n2win is None or n3win is None:
		n1win=n1;
		n2win=n2;
		n3win=n3;
	
	nov1=(1-r1)*n1win; nov1=int(nov1); # non-overlapping size 1
	nov2=(1-r2)*n2win; nov2=int(nov2); # non-overlapping size 2
	nov3=(1-r3)*n3win; nov3=int(nov3); # non-overlapping size 3
	ov1=r1*n1win;      ov1=int(ov1);   # overlapping size 1
	ov2=r2*n2win;      ov2=int(ov2);  # overlapping size 2
	ov3=r3*n3win;      ov3=int(ov3);  # overlapping size 3

	n1pad=n1win;        # padding size 1
	nw1=1;
	while n1pad<n1:
		n1pad=n1pad+nov1;nw1=nw1+1;


	n2pad=n2win;        # padding size 2
	nw2=1;
	while n2pad<n2:
		n2pad=n2pad+nov2;nw2=nw2+1;

	n3pad=n3win;        # padding size 3
	nw3=1;
	while n3pad<n3:
		n3pad=n3pad+nov3;nw3=nw3+1;


	D1=np.zeros([n1pad,n2pad, n3pad]);D1[0:n1,0:n2,0:n3]=din;#copy din into padded D1
	mask_tmp=mask;mask=np.zeros([n1pad,n2pad,n3pad]);mask[0:n1,0:n2,0:n3]=mask_tmp;
	D2=np.zeros([n1pad,n2pad,n3pad]);

	for iw3 in range(0,nw3):
		for iw2 in range(0,nw2):
			for iw1 in range(0,nw1):
				s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;
				dtmp=D1[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win];
				
				# uncomment this line for checking the correctness (checked 100%  correct)
				param['mask']=mask[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win];
				dtmp = oper(dtmp,param);
				# only valid for space-independent param
				# for reconstruction, with mask, param needs to be changed
				
				dtmp=win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
				D2[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win]=D2[s1:s1+n1win,s2:s2+n2win,s3:s3+n3win]+dtmp;
	dout=D2[0:n1,0:n2,0:n3];
	return dout
	

