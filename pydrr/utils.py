def scale(D,N=2,dscale=1.0):
	"""
	scale: Scale the data up to the Nth dimension = sfscale axis=N
	IN   D:   	intput data
	     N:      number of dimension for scaling
	             default: N=2
		dscale:  Scale by this factor
	    (does not include the rscale and pclip functions (not convenient actually))
 
	OUT   D1:  	output data
	
	Copyright (C) 2015 The University of Texas at Austin
	Copyright (C) 2015 Yangkang Chen
	Modified by Yangkang Chen on Jan, 2020
 
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
 
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	"""
	import numpy as np

	if D.ndim==2:	#for 2D problems
		D=np.expand_dims(D, axis=2)
		
	if D.ndim==1:	#for 1D problems
		D=np.expand_dims(D, axis=1)
		D=np.expand_dims(D, axis=2)
		
	[n1,n2,n3]=D.shape;
	
	D1=D;
	
	if N==1:
		for i3 in range(0,n3):
			for i2 in range(0,n2):
				D1[:,i2,i3]=D1[:,i2,i3]/np.max(np.abs(D1[:,i2,i3]));
	elif N==2:
		for i3 in range(0,n3):
			D1[:,:,i3]=D1[:,:,i3]/np.max(np.abs(D1[:,:,i3]));
	elif N==3:
		D1=D1/np.max(np.abs(D1));
	elif N==0:
		D1=D1*dscale;
	else:
		print("Invalid argument value N.");

	
	
	D1=np.squeeze(D1);
	
	return D1