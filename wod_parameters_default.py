import numpy as np
import healpy as hp
pi = np.pi

m_lim = 19.8
nside = 256 
nside_tag = '_nside' + str(nside)
nside_jk = 8 
N_pix = hp.nside2npix(nside)

ang_min = 0.45
ang_max = 10. 
ang_res = 0.3

cat_list = {'p':True, 's':True}
map_list = {'star':True, 'gal':True, 'odds':True}
correction_list = ['od', 'star'] 

cat = {} 
if(cat_list['s'] == True): cat['s'] = {'val':{}}
if(cat_list['p'] == True): cat['p'] = {'val':{}}

cat['s']['pwd'] = '/Users/pmarti/Dropbox/Tesi/photoz/2slaq_DR7_LRG_webcat/BPZ/'
cat['p']['pwd'] = '/Users/pmarti/Dropbox/Tesi/photoz/MegaZ_DR7/BPZ/'

cat['s']['f'] = '2slaq_DR7_LRG_webcat_hdr_nonobserved.bpz' 
cat['p']['f'] = 'MegaZ_DR7_i19.8.bpz' 

mask_file = '/Users/pmarti/Dropbox/Tesi/correlation/MegaZ_DR7/map/MegaZ_DR7_nside256_anne_mask_map.fits'

star_map_file = '/Users/pmarti/Dropbox/Tesi/correlation/MegaZ_DR7/map/stars_r9.fits'

cat['s']['col'] = {'zp': 1, 'zs': 10, 'od':5}
cat['p']['col'] = {'zp': 1, 'od':5, 'ra':11, 'dec':12}

#zbin = [[0.45,0.5],[0.5,0.55],[0.55,0.6],[0.6,0.65]]
zbin = [[0.5,0.6]]
n_zbin = len(zbin)

#od_cut = np.arange(0,0.6,0.1)
od_cut = np.array([0.]) 
n_od = len(od_cut)
od_range = [0.4,1.]
