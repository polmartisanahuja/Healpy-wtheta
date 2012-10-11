import os
import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
pi = np.pi
from parameters12 import *

#MAIN..........................................................
print "******************************************"
print "*                                        *"
print "*       Healpix Map12 by Polstein        *" 
print "*                                        *"
print "******************************************"

#Creating zbin folders
nz_dir_name = './Nz/' + dir_name
plot_dir_name = './plot/' + dir_name
if not os.path.exists(nz_dir_name):
	os.makedirs(nz_dir_name)
if not os.path.exists(plot_dir_name):
	os.makedirs(plot_dir_name)

#Printing input parameters
print "\nPhotoz bins:"
print "\tLow zbin: from %.2f to %.2f" % (Zbin[0][0], Zbin[0][1])
print "\tHigh zbin: from  %.2f to %.2f" % (Zbin[1][0], Zbin[1][1])
print "\nOdds cut list:"
for od_cut in Od_cut: print "\tod > %.2f" % od_cut

from map_functions import *

#Loading catalog
print "\nLoading " + theo_cat_name + " catalog..." 
theo_cat_origin = np.loadtxt(theo_cat_file, usecols = theo_columns, unpack = True)
print "\nLoading " + cat_name + " catalog..." 
cat_origin = np.loadtxt(cat_file, usecols = columns, unpack = True)

#Creating mask
if(os.path.isfile(mask_file)): mask_map = hp.read_map(mask_file)
else:
	cat_mask = load_cat(cat_origin, zbin = False, od_cut = False)
	mask_map = mask_map(cat_mask)
	hp.write_map(mask_file, mask_map)

print "\nComputing maps and Nz:" 
leg = True 
for zbin in Zbin:
	c = 0
	for od_cut in Od_cut:

		#Tags and names
		zbin_tag = "_zbin%.2f-%.2f" % (zbin[0], zbin[1])	
		od_tag = "_od%.2f" %  od_cut
		name = cat_name + nside_tag  + zbin_tag + od_tag

		print "\n\t*******N Map " + zbin_tag + od_tag + "*******"
	
		#Load catalog
		cat = load_cat(cat_origin, zbin, od_cut) 
		
		#n map
		n_map = dens_map(cat)
		n_map_file = './map/' + name + "_n_map.fits"
		hp.write_map(n_map_file, n_map)
				
		#Cuting theoretical catalog
		theo_cat = cut(theo_cat_origin, theo_zphot_col, "zp", zbin[0], zbin[1])
		theo_cat = cut(theo_cat, theo_odds_col, "od", od_cut, 1.)

		#Histogram
		Nz, z = np.histogram(theo_cat[theo_zspec_col], range = (bin_range[0], bin_range[1]), bins = nbins)
		z = z[:-1] + bin_size / 2
		
		#Savefile
		nz_file = nz_dir_name + '/' + theo_cat_name + zbin_tag + od_tag + '_Nz.txt'  
		np.savetxt(nz_file, np.array([z,Nz]).T, fmt = '%.3f %d')

		#Plot
		plt.plot(z, Nz, color = colors[c], label = 'odds > %.2f' % od_cut, linewidth = 2)
		c += 1

		if(od_cut == 0.):
			#Odds map
			print "\n\t*******Od Map" + zbin_tag + "*******"
			od_map = odd_map(cat, mask_map)
			od_map_file = './map/' + name + "_od_map.fits"
			hp.write_map(od_map_file, od_map)

	if(leg): plt.legend()
	leg = False

plt.title(theo_cat_name + ' zbin' + z_tag[0] + '_' + z_tag[1])
plt.ylabel('Counts')
plt.xlabel('z(spec)')
#plt.show()
plt.savefig(plot_dir_name + '/nz.png')
