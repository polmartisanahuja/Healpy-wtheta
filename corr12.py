import os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import healpy as hp
import math
import time
pi = np.pi
from parameters12 import *

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
		
		corr_ij['n12'][k] += n_map[1][i] * n_map_masked[0][disc].sum()
		corr_ij['od12'][k] += od_map[1][i] * od_map_masked[0][disc].sum()
		corr_ij['od1'][k] += od_map[0][i] * od_map_masked[0][disc].sum()
		corr_ij['od2'][k] += od_map[1][i] * od_map_masked[1][disc].sum()
		corr_ij['nod1'][k] += n_map[0][i] * od_map_masked[0][disc].sum()
		corr_ij['nod2'][k] += n_map[1][i] * od_map_masked[1][disc].sum()
		n_ij[k] += mask_map[disc].sum()
	
	#Take non-zero elements
	for l in corr_list[:-1]: corr_ij[l] = corr_ij[l][np.nonzero(n_ij)]
	n_ij = n_ij[np.nonzero(n_ij)]	
	K = len(corr_ij['n12'])

	#Mean correlations
	mean = {}
	for l in corr_list[:-1]: mean[l] = corr_ij[l].sum() / n_ij.sum()
	mean['n12_corrected'] = mean['n12'] - (mean['nod1'] * mean['nod2'] * mean['od12']) / (mean['od1'] * mean['od2']) 

	#Jackknife correlations
	kcorr_ij = {} 
	for l in corr_list[:-1]:
		kcorr_ij[l] = np.empty(K)
		for k in range(K):
			kcorr_ij[l][k] = (corr_ij[l].sum() - corr_ij[l][k]) / (n_ij.sum() - n_ij[k])
	kcorr_ij['n12_corrected'] = np.empty(K)

	for k in range(K): kcorr_ij['n12_corrected'][k] = kcorr_ij['n12'][k] - (kcorr_ij['nod1'][k] * kcorr_ij['nod2'][k] * kcorr_ij['od12'][k]) / (kcorr_ij['od1'][k] * kcorr_ij['od2'][k])
		
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
print "*       Correlations12 by P.Marti        *" 
print "*                                        *"
print "******************************************"

#Printing parameters
print "\nAng max = ", ang_max
print "Ang res = ", ang_res
print "Nside = ", nside
print "Nside jk = ", nside_jk
print "zbin1 = ",  z_tag[0] 
print "zbin1 = ",  z_tag[1]
print "Odds cut list = ", Od_cut 
print "Mask file = ", mask_file

#Create output directory
#dir_name = './corr/zbin%.2f-%.2f_%.2f-%.2f' % (Zbin[0][0], Zbin[0][1], Zbin[1][0], Zbin[1][1])
#dir_name = './corr/zbin' + z_tag[0] + '_' + z_tag[1]
corr_dir_name = './corr/' + dir_name
if not os.path.exists(corr_dir_name):
	os.makedirs(corr_dir_name)

#Loading mask
mask_map = hp.read_map(mask_file)

#Getting pixel index in the mask
mask = np.nonzero(mask_map)[0]

#Computing correlations in parallel
ang = np.arange(1.5 * ang_res, ang_max, ang_res)

n_map = np.empty([2,N])
n_map_masked = np.empty([2,N])
od_map = np.empty([2,N])
od_map_masked = np.empty([2,N])
#zbin_tag = np.empty([2], dtype = np.object_)

for od_cut in Od_cut:
	
	print "\n***********Odds cut = %.2f***********" % od_cut
	od_tag = "_od%.2f" %  od_cut
	
	for i in range(2):
		#zbin_tag[i] = "_zbin%.2f-%.2f" % (Zbin[i][0], Zbin[i][1])
		zbin_tag = '_zbin' + z_tag[i]	
		name = cat_name + nside_tag  + zbin_tag + od_tag
		n_map_file = './map/' + name + '_n_map.fits'
		n_map[i] = hp.read_map(n_map_file)
		mean_n_map = (n_map[i] * mask_map).sum() / mask_map.sum()
		print "mean_n_map%d = %s" %  (i + 1, mean_n_map)
		n_map[i] = (n_map[i] / mean_n_map) - 1. 
		n_map_masked[i] = mask_map * n_map[i]

		if(od_cut == 0.): 
			od_map_file = './map/' + name + '_od_map.fits'
			od_map[i] = hp.read_map(od_map_file)
			mean_od_map = (od_map[i] * mask_map).sum() / mask_map.sum()
			print "mean_od_map%d = %s" %  (i + 1, mean_od_map)
			od_map[i] = (od_map[i] / mean_od_map) - 1. 
			od_map_masked[i] = mask_map * od_map[i]
	
	#correlations(ang[0])
	#import sys
	#sys.exit()

	po = Pool() 
	T0 = time.time()
	res = po.map(correlations, ang)
	T = time.time() - T0
	print "\tTotal time = %1.1f sec" % T 

	for l in corr_list:
		corr = []
		err_corr = []
		for val in res: 
			corr = np.append(corr, val[0][l])
			err_corr = np.append(err_corr, val[1][l])

		#Writting correlations 
		corr_file = corr_dir_name + '/' + cat_name + nside_tag + nside_jk_tag + '_zbin' + z_tag[0] + '_' + z_tag[1] + od_tag + '_' + l +".txt" 
		np.savetxt(corr_file, np.array([ang, corr, err_corr]).T)
