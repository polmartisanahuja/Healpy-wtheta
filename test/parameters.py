from enviromental import *
import numpy as np
cat_folder = "./" 
cat_name = "2slaq_DR7_LRG_webcat_hdr"
cat_file = cat_folder + cat_name + ".bpz"
columns = (1,5,11,12)
odds_col = 1
zphot_col = 0
ra_col = 2
dec_col = 3

ang_max = 10.  #Maximum angle of the correlations
ang_res = 0.3 #Resolution of the angle bins in the correlations

nside = 256 #Parameter that defines the resolution of the healpix map (N_pix = 12 * nside^2)
nside_jk = 8 #Resolution of the jackknife pixels

z_cut_low = 0.5 #Photo-z bin edges
z_cut_high = 0.6

odds_cut_low = 0.65 #Odds cut 
odds_cut_high = 1.0 
od = "%.2f" % odds_cut_low
odds_cut_list = [0.00, 0.65]

corr_list = np.array(['n_auto', 'od_auto', 'cross', 'n_auto_corrected'])

n_map_file = cat_name + "_nside" + str(nside) + "_n_map_od" + od + ".fits"
od_map_file = cat_name + "_nside" + str(nside) + "_odds_map.fits"
mask_file = cat_name + "_nside" + str(nside) + "_mask.fits"

cord_map = cat_name + "_ra&dec.txt"

#map1_file = n_map_file
#map2_file = n_map_file

#corr_type = "auto_gal_"
