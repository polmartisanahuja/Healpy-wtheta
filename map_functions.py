import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
pi = np.pi

from parameters12 import *

#Printing INPUT PARAMETERS
pixresol = hp.nside2resol(nside) #Radians
N_pix = hp.nside2npix(nside)
S_pix = hp.nside2pixarea(nside, degrees = True)

print "\nHealpix Parameters:"
print "\tNSIDE = %d" % nside
print "\t# total pixels = %d" % N_pix
print "\tPixel resolution = %3.3f (deg)" % ((180. / pi) * pixresol)
print "\tPixel area = %3.3f (sq.deg.)" % S_pix 
print "\tPixel size = %3.3f (Mpc) (LCDM cosmology)" % (pixresol * 1254.5)

def cut(cat, col, name,  cut_low, cut_high):

	x = cat[col]
	N = len(x) #N: Total number of galaxies
	rm_x_pos = []
	for i in range(len(x)):
		if ((cut_low <= x[i]) & (x[i] <= cut_high)): rm_x_pos.append(i)
	cat = np.take(cat, rm_x_pos, 1)

	#Print proportions.................................
	n = len(cat[col]) #n: Galaxies after cut
	print "\n\t\t" + name + " cut:"
	print "\t\t# gal = %d" % N
	print "\t\t# gal. after $%2.2f<%s<%2.2f$ cut = %d %1.1f%%" % (cut_low, name, cut_high, n , 100 * float(n)/float(N) )
	return cat

def load_cat(cat, zbin, od_cut):

	#Loading catalog
	#cat = np.loadtxt(cat_file, usecols = columns, unpack = True)

	#zbin cut
	if(zbin): cat = cut(cat, zphot_col, "zp", zbin[0], zbin[1])

	#Odds cut
	if(od_cut): cat = cut(cat, odds_col, "od", od_cut, 1.)

	#Conversion deg to rad
	cat[ra_col] *= (pi / 180)
	cat[dec_col] *= (pi / 180)

	#Changing from Equatorial coordinates to HEALPix coordinates
	phi = cat[ra_col] 
	theta = (pi / 2) - cat[dec_col] 

	#Coverting angle coord to pixel index
	pix = hp.ang2pix(nside, theta, phi)

	#Creating new cat array with pix index instead of ra,dec
	cat = np.append(cat[:-2], [pix], axis = 0)	

	return cat

def dens_map(cat):

	n_map = np.zeros(N_pix)
	pix = cat[-1]
	for i in pix: n_map[i] += 1
	
	return n_map

def mask_map(cat):
	"""Returns an histogram of pixels within the MASK.
	Tricky MASK!!: By now, this is a provisional mask that
	that removes N>5 contiguos pixels with zeros content."""

	#Loading n_map
	n_map = dens_map(cat)

	pix_index = np.nonzero(n_map)[0]
	N = len(pix_index)

	#Find zeros_threshold
	print "Consecutive zeros histogram:"
	D = []
	for i in range(1, N): 
		sys.stdout.write("\r\tComplete: %1.1f%%" % (100 * float(i) / float(N)))
		sys.stdout.flush()
		D = np.append(D, pix_index[i] - pix_index[i - 1])

	n_bins = 30
	plt.hist(D, bins = n_bins, range = (2, n_bins + 2))
	plt.xlabel("consecutive zeros")
	plt.show()

	zeros_trsh = input("\nzeros_trsh in the MASK:")

	#Re-add zero content pixels due to statistical fluctiations
	print "Re-add statistical zero chains:"
	zeros_index = []
	zeros_pix_index = []
	for i in range(1,N):
		sys.stdout.write("\r\tComplete: %1.1f%%" % (100 * float(i) / float(N)))
		sys.stdout.flush()
		d = pix_index[i] - pix_index[i - 1]
		if ((d < zeros_trsh) & (d > 1)): 
			zeros_index = (np.append(zeros_index, i * np.ones(d - 1))).astype(int)
			zeros_pix_index = (np.append(zeros_pix_index, np.arange(pix_index[i - 1] + 1,  pix_index[i], 1,))).astype(int)
	pix_index = np.insert(pix_index, zeros_index, zeros_pix_index )

	#Create the mask throught pix index
	mask_map = np.zeros(hp.nside2npix(nside))
	for i in pix_index:
		mask_map[i] = 1

	return mask_map 

def odd_map(cat, mask_map):
	
	#Loading mask
	#mask_map = hp.read_map(mask_file)

	#Loading n_map
	n_map = dens_map(cat)
	
	#Creating od_map
	odds_map = np.zeros(N_pix)
	odds_mean = cat[odds_col].sum() / len(cat[odds_col]) #Odds mean of the catalog
	print "\t\tOdds cat mean  = %4.4f" % odds_mean 	

	pix = cat[-1]
	for i in range(len(pix)): odds_map[pix[i]] += cat[odds_col][i]  
	
	n_noempty = 0 #Number of non-empty pixels in the mask
	mask = np.nonzero(mask_map)[0] #Index of the pixels in the mask
	n_mask = len(mask)
	for i in mask: 
		if(n_map[i] != 0): 
			odds_map[i] /= n_map[i] 
			n_noempty += 1

	print "\t\tFraction of empty pixels = %4.4f%%" % (100*float(n_mask - n_noempty) / float(n_mask)) 	

	#Correcting od_map empty pixels (1deg averaged)
	n_after = 0
	odds_map_corrected = np.copy(odds_map)
	for i in mask: 
		if(n_map[i] == 0): 
			i_vec = hp.pix2vec(nside, i)
			nest = hp.query_disc(nside, i_vec, 1.0*(np.pi/ 180))
			odds = 0
			n = 0
			for j in nest: 
				odds += odds_map[j]
				if (odds_map[j] > 0.00000001): n += 1
			odds /= float(n)
			odds_map_corrected[i] = odds
			if (odds < 0.00000001): n_after += 1 

	print "\t\tOdds map mean = %4.4f" % (odds_map.sum() / n_mask) 	
	print "\t\tCorrected odds map mean = %4.4f" % (odds_map_corrected.sum() / n_mask) 	

	return odds_map_corrected
