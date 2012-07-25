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
def plot_corr():
	"""Generates a plot of the angular correlations"""

	plt.close()
	plt.errorbar(ang, corr, err_corr)
	plt.axhline( 0, linewidth = 1, color = 'red')
	plt.xlabel("angle (deg)")
	plt.ylabel("<map1 * map2>")
	plt.savefig(corr_plot)

def plot_time():
	""" Generates a plot of the computation time for each angle bin"""
	plt.close()
	plt.plot(ang, t, 'o-')
	plt.xlabel("angle (deg)")
	plt.ylabel("Time (seconds)")
	plt.savefig(time_plot)

def correlations(ang):
	"""Computes correlations of delta in an angle bin of ang_bin_res resolution. Input angle in radians"""
	
	n_ij = np.zeros(N_jk)
	sum_ij = np.zeros(N_jk)
	
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

		for j in disc:
			sum_ij[k] += map1[i] * map2[j] * mask_map[j]
			n_ij[k] += mask_map[j] 
	
	#Take non-zero elements
	sum_ij = sum_ij[np.nonzero(sum_ij)]	
	n_ij = n_ij[np.nonzero(n_ij)]	
	K = len(sum_ij)

	#Jackknife correlations
	corr = np.empty(K)
	for k in range(K): corr[k] = (sum_ij.sum() - sum_ij[k]) / (n_ij.sum() - n_ij[k])
	
	#Mean correlations
	mean_corr =  sum_ij.sum() / n_ij.sum()
	
	#Error correlation
	err_corr = 0
	for k in range(K): err_corr += (corr[k] - mean_corr) * (corr[k] - mean_corr)
	err_corr = math.sqrt(((float(K) - 1.) / float(K)) * err_corr)

	t = time.clock() - t0

	print "\tAng, Corr, Err_corr = %4.4f, %7.7f, %7.7f" % (ang, mean_corr, err_corr)
	return mean_corr, err_corr, t

#MAIN..........................................................

print "******************************************"
print "*                                        *"
print "*         Correlations by P.Marti        *" 
print "*                                        *"
print "******************************************"

#Loading maps 
map1 = hp.read_map(map1_file)
map2 = hp.read_map(map2_file)

#Loading mask
mask_map = hp.read_map(mask_file)

#Getting pixel index in the mask
mask = np.nonzero(mask_map)[0]

mean_map1 = (map1 * mask_map).sum() / mask_map.sum()
mean_map2 = (map2 * mask_map).sum() / mask_map.sum()

print "mean_map1=", mean_map1
print "mean_map2=", mean_map2

map1 = (map1 / mean_map1) - 1. 
map2 = (map2 / mean_map2) - 1. 

#Computing correlations in parallel
ang = np.arange(1.5 * ang_res, ang_max, ang_res)

po = Pool() 
T0 = time.time()
res = po.map(correlations, ang)
T = time.time() - T0
print "\tTotal time = %1.1f sec" % T 

corr = []
err_corr = []
t = []
for val in res: 
	corr = np.append(corr, val[0])
	err_corr = np.append(err_corr, val[1])
	t = np.append(t, val[2])

#Plots
plot_corr()
plot_time()

#Writting correlations 
np.savetxt(corr_file, np.array([ang, corr, err_corr]).T)
