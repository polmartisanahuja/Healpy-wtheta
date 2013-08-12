import numpy as np
import healpy as hp
pi = np.pi

def zbintag(bin): return 'zbin' + str(bin[0]) + '-' + str(bin[1])

def odtag(od): return 'od' + str(od)

def angtag(ang_min, ang_max): return '_ang' + str(ang_min) + '-' + str(ang_max)

def correctiontag(l): 
	if(len(l) == 0): return ''
	else:
		s = '_' 
		for i in l: s += i + '_'
		return s + 'corrected' 

def cut(zbin_cut, od_cut, zp, od):
	
	mask_od = od >= od_cut 
	mask_bin = (zp >= zbin_cut[0]) & (zp <= zbin_cut[1])
	return np.logical_and(mask_od, mask_bin)

def idmaskmap(mask_map):

	N = len(mask_map)
	idmask_map = np.zeros(N, dtype = 'int')
	j = 0
	for i in range(N):
		if(mask_map[i] == 0): idmask_map[i] = -1 
		else: 
			idmask_map[i] = j 
			j += 1
	return idmask_map

def odrescal(od_array, od_range):

	od_array = (od_array - od_range[0]) / (od_range[1] - od_range[0])
	mask = od_array < 0
	od_array[mask] = 0
	return od_array

def paths(tree, cur=()):
	if type(tree) != type({}):
		yield cur
	else:
		for n, s in tree.items():
			for path in paths(s, cur+(n,)):
				yield path

def delta(map):
	
	mean = map.sum() / len(map)
	return (map / mean) - 1.

def coord2pix(ra, dec, nside):

	#Conversion deg to rad
	ra *= (pi / 180)
	dec *= (pi / 180)

	#Changing from Equatorial coordinates to HEALPix coordinates
	phi = ra 
	theta = (pi / 2) - dec 

	#Coverting angle coord to pixel index
	return hp.ang2pix(nside, theta, phi)

def get_dict(c,l):
	
	for i in l: c = c[i] 
	return c

def corrtag(l):
	
	tag = ''
	for i in l: 
		tag += "_" + i	
	return tag

def correction_eq(n, ns1, ns2, s12, s1, s2): return n - (ns1 * ns2 * s12) / (s1 * s2)
