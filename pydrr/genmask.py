import numpy as np
def genmask(u, r, type, seed):
	"""
	GENMASK:Generate Random Sampling Mask
	
	INPUT
	u: 		image
	r: 		data KNOWN ratio
	type: 	data lose type
	   		'r': random lose rows
	   		'c': random lose columns
	   		'p': random lose pixel
	seed: 	seed of random number generator
	
	OUTPUT
	mask: 	sampling mask
	
	  Copyright (C) 2014 The University of Texas at Austin
	  Copyright (C) 2014 Yangkang Chen
	  Ported 
	
	  This program is free software: you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation, either version 3 of the License, or
	  any later version.
	
	  This program is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	  GNU General Public License for more details:
	  http://www.gnu.org/licenses/
	"""
	
	m=u.shape[0];
	n=u.shape[1];
	

	mask = np.zeros([m,n]);
	
	if type=='r':
		row = rperm(m,seed);
		k = np.fix(r*m);k=int(k);
		row = row[0:k-1];
		mask[row,:] = 1;
		
	elif type=='c':
		column = rperm(n,seed);
		k = np.fix(r*n);k=int(k);
		column = column[0:k-1];
		mask[:, column] = 1;
	elif type=='p':
		pix = rperm(m*n,seed);
		k = np.fix(r*m*n);k=int(k);
		pix = pix[0:k-1];
		mask[pix]= 1;
	else:
		print("mask type not found");
		
	return mask

def rperm(n,seed):
	"""
	RPERM: Random permutation of my version.
	
	RPERM(n) is a random permutation of the integers from 1 to n.
	For example, RANDPERM(6) might be [2 4 5 6 1 3].
	
	  Copyright (C) 2014 The University of Texas at Austin
	  Copyright (C) 2014 Yangkang Chen
	
	  This program is free software: you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation, either version 3 of the License, or
	  any later version.
	
	  This program is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	  GNU General Public License for more details:
	  http://www.gnu.org/licenses/
	"""
	
	np.random.seed(seed);
	p = np.argsort(np.random.rand(n));

	return p
