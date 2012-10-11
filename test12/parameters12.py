from enviromental import *
import numpy as np
import healpy as hp

#Input catalog
cat_folder = "./" 
cat_name = "2slaq_DR7_LRG_webcat_hdr" 
cat_file = cat_folder + cat_name + ".bpz"
columns = (1,5,11,12)
odds_col = 1
zphot_col = 0
ra_col = 2
dec_col = 3

#Input theoretical catalog
theo_cat_folder = "./" 
theo_cat_name = "2slaq_DR7_LRG_webcat_hdr"
theo_cat_file = theo_cat_folder + theo_cat_name + ".bpz" 
theo_columns = (1,5,10)
theo_odds_col = 1
theo_zphot_col = 0
theo_zspec_col = 2

#Angle resolution and range
ang_max = 10.  #Maximum angle of the correlations
ang_res = 0.3 #Resolution of the angle bins in the correlations

#Map resolution
nside = 256 
N = hp.nside2npix(nside)
nside_tag = "_nside" + str(nside)

#Jacknife resolution
nside_jk = 8 
nside_jk_tag = "_jk" + str(nside_jk)
N_jk = hp.nside2npix(nside_jk)

#Photo-z bin edges
Zbin = [[0.50, 0.60], [0.50, 0.60]]
z_tag = np.empty([2], dtype = np.object_)
for i in range(2): z_tag[i] = "%.2f-%.2f" % (Zbin[i][0], Zbin[i][1])

#Odds cut
Od_cut = [0.00, 0.55, 0.65]

corr_list = np.array(['n12', 'od12', 'od1','od2', 'nod1','nod2', 'n12_corrected'])

#File names
mask_file = "./map/" + cat_name + nside_tag + "_mask_map.fits"
dir_name = 'zbin' + z_tag[0] + '_' + z_tag[1]

#Histogram parameters
bin_size = 0.01
bin_range = [0.2,0.9]
nbins = (bin_range[1] - bin_range[0]) / bin_size

#Plot on/off
theo_curves = False
colors = np.array(['g', 'b', 'r', 'c', 'm', 'y', 'k'])
pn = True
pod = True
pnod = True
log_scale = False
y_lim = [-0.05, 0.14]
#y_lim = False 
#x_lim = [0., 10.]
x_lim = False 
