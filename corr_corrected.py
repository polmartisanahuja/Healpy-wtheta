import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import healpy as hp
import math
import time

#Parameters...................................................
from parameters import *
N_jk = hp.nside2npix(nside_jk)
hp_version = hp.__version__
pi = np.pi

#Functions....................................................
def correlations(ang):
	"""Computes correlations of delta in an angle bin of ang_bin_res resolution. Input angle in radians"""
	
	n_ij = np.zeros(N_jk)
	corr_ij = {}
	for l in corr_list[:-1]: corr_ij[l] = np.zeros(N_jk) 
	
	ang_low = (pi / 180) * (ang - (ang_res / 2))
	ang_high = (pi / 180) * (ang + (ang_res / 2))

	t0 = time.clock()
	for i in mask:
		#Get pixel in a disc of angle ang and width ang_res arround i
		i_vec = hp.pix2vec(nside, i) #Convert pixel index to unitary vector
		disc_low = hp.query_disc(nside, i_vec, ang_low, inclusive = False) #Pixels inside smaller disc 
		disc_high = hp.query_disc(nside, i_vec, ang_high, inclusive = False) #Pixels inside bigger disc
		disc = np.setdiff1d(disc_high, disc_low) #Takes non-overlaped pixels 
		
		k = hp.vec2pix(nside_jk, i_vec[0], i_vec[1], i_vec[2]) #Looking for which jackknife pixel corresponds to
		
		corr_ij['n_auto'][k] += n_map[i] * n_map_masked[disc].sum()
		corr_ij['od_auto'][k] += od_map[i] * od_map_masked[disc].sum()
		corr_ij['cross'][k] += n_map[i] * od_map_masked[disc].sum()
		n_ij[k] += mask_map[disc].sum()
	
	#Take non-zero elements
	for l in corr_list[:-1]: corr_ij[l] = corr_ij[l][np.nonzero(n_ij)]
	n_ij = n_ij[np.nonzero(n_ij)]	
	K = len(corr_ij['n_auto'])

	#Mean correlations
	mean = {}
	for l in corr_list[:-1]: mean[l] = corr_ij[l].sum() / n_ij.sum()
	mean['n_auto_corrected'] = mean['n_auto'] - mean['cross'] * mean['cross'] / mean['od_auto'] 

	#Jackknife correlations
	kcorr_ij = {} 
	for l in corr_list[:-1]:
		kcorr_ij[l] = np.empty(K)
		for k in range(K):
			kcorr_ij[l][k] = (corr_ij[l].sum() - corr_ij[l][k]) / (n_ij.sum() - n_ij[k])
	kcorr_ij['n_auto_corrected'] = np.empty(K)
	

	for k in range(K): kcorr_ij['n_auto_corrected'][k] = kcorr_ij['n_auto'][k] - kcorr_ij['cross'][k] * kcorr_ij['cross'][k] / kcorr_ij['od_auto'][k] 
		
	#Error correlation
	err_corr = {}
	print "\n#Theta = %4.4f" % ang
	for l in corr_list:
		err_corr[l] = 0
		for k in range(K): err_corr[l] += (kcorr_ij[l][k] - mean[l]) * (kcorr_ij[l][k] - mean[l]) 
		err_corr[l] = math.sqrt(((float(K) - 1.) / float(K)) * err_corr[l])
		print "%s: %7.7f, %7.7f" % (l, mean[l], err_corr[l])

	t = time.clock() - t0

	return mean, err_corr, t

#MAIN..........................................................

print "******************************************"
print "*                                        *"
print "*         Correlations by P.Marti        *" 
print "*                                        *"
print "******************************************"

#Loading maps 
n_map = hp.read_map(n_map_file)
od_map = hp.read_map(od_map_file)

#Loading mask
mask_map = hp.read_map(mask_file)

#Getting pixel index in the mask
mask = np.nonzero(mask_map)[0]

mean_n_map = (n_map * mask_map).sum() / mask_map.sum()
mean_od_map = (od_map * mask_map).sum() / mask_map.sum()

print "mean_n_map=", mean_n_map
print "mean_od_map =", mean_od_map

n_map = (n_map / mean_n_map) - 1. 
od_map = (od_map / mean_od_map) - 1. 

n_map_masked = mask_map * n_map
od_map_masked = mask_map * od_map

#Computing correlations in parallel
ang = np.arange(1.5 * ang_res, ang_max, ang_res)

po = Pool() 
T0 = time.time()
res = po.map(correlations, ang)
T = time.time() - T0
print "\tTotal time = %1.1f sec" % T 

for l in corr_list:
	corr = []
	err_corr = []
	t = []
	for val in res: 
		corr = np.append(corr, val[0][l])
		err_corr = np.append(err_corr, val[1][l])
		t = np.append(t, val[2])

	#Writting correlations 
	corr_file = "corr/" + l + "_od" + od + ".txt" 
	np.savetxt(corr_file, np.array([ang, corr, err_corr]).T)
