import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from math import *
pi = np.pi

#Parameters..................................................
from parameters import *

#Functions.....................................................
def cut(cat, col, name,  cut_low, cut_high):

	x = cat[col]
	N = len(x) #N: Total number of galaxies
	rm_x_pos = []
	for i in range(len(x)):
		if ((cut_low <= x[i]) & (x[i] <= cut_high)): rm_x_pos.append(i)
	cat = np.take(cat, rm_x_pos, 1)

	#Print proportions.................................
	n = len(cat[col]) #n: Galaxies after cut
	print "\n" + name + " cut:"
	print "\tNum. gal = %d" % N
	print "\tNum. gal. after $%2.2f<%s<%2.2f$ cut = %d %1.1f%%\n" % (cut_low, name, cut_high, n , 100 * float(n)/float(N) )
	return cat

#MAIN..........................................................

print "******************************************"
print "*                                        *"
print "*         Healpix Map by P.Marti         *" 
print "*                                        *"
print "******************************************"

#Loading catalog
cat = np.loadtxt(cat_file, usecols = columns, unpack = True)

#Redshift cut
cat = cut(cat, zphot_col, "zphot", z_cut_low, z_cut_high)

#Odds cut
#cat = cut(cat, odds_col, "odds", odds_cut_low, odds_cut_high)

#Magnitude cut
#cat = cut(cat, ideV_col, "ideV",  m_cut_low, m_cut_high)

#Coordinates file
np.savetxt(cord_map, np.array([cat[ra_col], cat[dec_col]]).T, fmt = '%8.8f %8.8f')

#Conversion deg to rad
cat[ra_col] *= (pi / 180)
cat[dec_col] *= (pi / 180)

#Changing from Equatorial coordinates to HEALPix coordinates
phi = cat[ra_col] 
theta = (pi / 2) - cat[dec_col] 

#Printing INPUT PARAMETERS
pixresol = hp.nside2resol(nside) #Radians
N_pix = hp.nside2npix(nside)
S_pix = hp.nside2pixarea(nside, degrees = True)

print "\nINPUT PARAMETERS:"
print "\tNSIDE = %d" % nside
print "\t# total pixels = %d" % N_pix
print "\tPixel resolution = %3.3f (deg)" % (180/ pi) * pixresol
print "\tPixel area = %3.3f (sq.deg.)" % S_pix 
print "\tPixel size = %3.3f (Mpc) (LCDM cosmology)" % (pixresol * 1254.5)

#Coverting angle coord to pixel index
pix = hp.ang2pix(nside, theta, phi)

#Galaxy number map..................................... 
#n_map, _ = np.histogram(pix, bins = N_pix, range = (0, N_pix))
n_map = np.zeros(N_pix)
for i in pix: n_map[i] += 1
hp.write_map(n_map_file, n_map)

#odds map.............................................. 
odds_mean = cat[odds_col].sum() / len(cat[odds_col])
odds_map = np.zeros(N_pix)
for i in range(len(pix)): odds_map[pix[i]] += cat[odds_col][i]  

#Count non-empty pixels in the mask (n)
mask_map = hp.read_map(mask_file)
mask = np.nonzero(mask_map)[0]
n_mask = len(mask)

n_noempty = 0
for i in mask: 
	if(n_map[i] != 0): 
		odds_map[i] /= n_map[i] 
		n_noempty += 1
odds_map_mean = odds_map.sum() / n_noempty 

print "\tOdds mean = %4.4f" % odds_mean 	
print "\tOdds mean map = %4.4f" % (odds_map.sum() / n_mask) 	
print "\tFraction of empty pixels = %4.4f%%" % (100*float(n_mask - n_noempty) / float(n_mask)) 	
print "\tCorrected odds mean map (whole mean) = %4.4f" % odds_map_mean 	

#Plot odds histogram
#plt.hist(cat[odds_col], bins = 99, range = (0.01, 1))
#plt.xlim(xmin = 0.01)
#plt.show()

#Plot odds map histogram
#plt.hist(np.take(odds_map, mask), bins = 99, range = (0.01, 1))
#plt.xlim(xmin = 0.01)
#plt.show()

#Odds-map correction 8nest
n_after = 0
odds_map_corrected = np.copy(odds_map)
for i in mask: 
	if(n_map[i] == 0): 
		i_vec = hp.pix2vec(nside, i)
		nest = hp.query_disc(nside, i_vec, 0.5*(np.pi/ 180))
		odds = 0
		n = 0
		for j in nest: 
			odds += odds_map[j]
			if (odds_map[j] > 0.00000001): n += 1
		odds /= float(n)
		odds_map_corrected[i] = odds
		if (odds < 0.00000001): n_after += 1 
print "\tCorrected odds map = %4.4f" % (odds_map_corrected.sum() / n_mask) 	
print "\tFraction of empty pixels =  %4.4f%%" % (100 * float(n_after) / float(n_mask)) 	

#Odds-map correction 8nest
#for i in mask: 
#	if(n_map[i] == 0): 
#		odds_map[i] = odds_map_mean 
#print "\tCorrected odds map (whole mean) = %4.4f" % (odds_map.sum() / n_mask) 	

#Plot corrected odds map histogram
plt.hist(np.take(odds_map_corrected, mask), bins = 99, range = (0.01, 1))
plt.xlim(xmin = 0.01)
plt.show()

#hp.write_map(var_map_file, odds_map)
hp.write_map(var_map_file, odds_map_corrected)
