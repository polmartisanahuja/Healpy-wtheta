import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from multiprocessing import Pool
import math
import time
import tools as tls
import os
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from wod_functions import *
from wod_megaz_para import *

#Functions......................................

def map_generator():

	#Constructing maps...............................
	#if not os.path.exists('./map'):
	#	os.makedirs('./map')

	#Star map
	if(map_list['star'] == True): map['star'] = hp.read_map(star_map_file)[mpix2hpix]

	#Galaxy map
	if((map_list['gal'] == True) or (map_list['odds'] == True)):
		for bin in zbin:
			map['gal'][zbintag(bin)] = {}
			i_eff = 0
			for od in od_cut:
				gal_map = np.zeros(N_mpix)	
				mask = cut(bin, od, cat['p']['val']['zp'], cat['p']['val']['od'])		
				mpix = cat['p']['val']['mpix'][mask]
				for i in mpix: 
					if(i >= 0): gal_map[i] += 1
				map['gal'][zbintag(bin)][odtag(eff_cut[i_eff])] = gal_map 
				
				if(map_list['gal'] == True): hp.write_map(folder_out + 'map/map_gal' + nside_tag + '_' + zbintag(bin) + '_' + odtag(eff_cut[i_eff]) + '.fits', np.append(gal_map, 0)[hpix2mpix])
				i_eff += 1

	#Odds map
	if(map_list['odds'] == True):
		for bin in zbin:
			mask = cut(bin, 0., cat['p']['val']['zp'], cat['p']['val']['od'])	
			od_map = np.zeros(N_mpix)	
			n_map = map['gal'][zbintag(bin)][odtag(0.0)]
			mpix = cat['p']['val']['mpix'][mask]
			odds = cat['p']['val']['od'][mask]
			for i in range(len(mpix)): 
				if(mpix[i] >= 0): od_map[mpix[i]] += odds[i] 
			for i in range(N_mpix):
				if(n_map[i] != 0): od_map[i] /= n_map[i]

			od_map_corr = np.copy(od_map)
			for i in range(N_mpix):
				if(n_map[i] == 0):
					i_vec = hp.pix2vec(nside, mpix2hpix[i])
					nest = hpix2mpix[hp.query_disc(nside, i_vec, 1 * (np.pi/ 180))]
					nest = nest[nest >= 0]
					o = 0
					n = 0
					for j in nest:
						o += od_map[j]
						if (od_map[j] > 0.00000001): n += 1
					o /= float(n)
					od_map_corr[i] = o
			
			map['odds'][zbintag(bin)] = od_map_corr
			hp.write_map(folder_out + 'Map/od_map' + nside_tag + zbintag(bin) + '.fits', np.append(od_map_corr,0)[hpix2mpix])
	
def set_corrdic():
	
	corr = {}
	
	for l in map_list: corr[l] = {}
	if(map_list['gal'] == True): corr['gal']['gal'] = {zbintag(zbin[i]):{zbintag(zbin[j]):{odtag(od):np.zeros((n_ang, N_jkmpix)) for od in eff_cut} for j in range(i, n_zbin)} for i in range(n_zbin)}
	if((map_list['gal'] == True) and (map_list['odds'] == True)): corr['gal']['odds'] = {zbintag(zbin[i]):{odtag(od):np.zeros((n_ang, N_jkmpix)) for od in eff_cut} for i in range(n_zbin)}
	if(map_list['odds'] == True): corr['odds']['odds'] = {zbintag(zbin[i]):{zbintag(zbin[j]):np.zeros((n_ang, N_jkmpix)) for j in range(i, n_zbin)} for i in range(n_zbin)}
	if((map_list['gal'] == True) and (map_list['star'] == True)): corr['gal']['star'] = {zbintag(zbin[i]):{odtag(od):np.zeros((n_ang, N_jkmpix)) for od in eff_cut} for i in range(n_zbin)}
	if((map_list['odds'] == True) and (map_list['star'] == True)): corr['odds']['star'] = {zbintag(zbin[i]):np.zeros((n_ang, N_jkmpix)) for i in range(n_zbin)}
	if(map_list['star'] == True): corr['star']['star'] = np.zeros((n_ang, N_jkmpix))	
	
	return corr

def get_map(l):
	if((l[0] == 'gal') and (l[1] == 'gal')): return map[l[0]][l[2]][l[4]], map[l[1]][l[3]][l[4]]
	if((l[0] == 'gal') and (l[1] == 'odds')): return map[l[0]][l[2]][l[3]], map[l[1]][l[2]]
	if((l[0] == 'gal') and (l[1] == 'star')): return map[l[0]][l[2]][l[3]], map[l[1]]
	if((l[0] == 'odds') and (l[1] == 'odds')): return map[l[0]][l[2]], map[l[1]][l[3]]
	if((l[0] == 'odds') and (l[1] == 'star')): return map[l[0]][l[2]], map[l[1]]
	if((l[0] == 'star') and (l[1] == 'star')): return map[l[0]], map[l[1]]


def write_corr(val,l):

	if(len(l) == 2): corr[l[0]][l[1]] = val 
	if(len(l) == 3): corr[l[0]][l[1]][l[2]] = val 
	if(len(l) == 4): corr[l[0]][l[1]][l[2]][l[3]] = val 
	if(len(l) == 5): corr[l[0]][l[1]][l[2]][l[3]][l[4]] = val  


def correction(l, select):
	
	n = corr['gal']['gal'][l[2]][l[3]][l[4]]

	if(select == 'od'):
		ns1 = corr['gal']['odds'][l[2]][l[4]]
		ns2 = corr['gal']['odds'][l[3]][l[4]]
		s12 = corr['odds']['odds'][l[2]][l[3]] 
		s1 = corr['odds']['odds'][l[2]][l[2]] 
		s2 = corr['odds']['odds'][l[3]][l[3]] 

	if(select == 'star'):
		ns1 = corr['gal']['star'][l[2]][l[4]]
		ns2 = corr['gal']['star'][l[3]][l[4]]
		s12 = corr['star']['star'] 
		s1 = corr['star']['star'] 
		s2 = corr['star']['star'] 
		
	write_corr(correction_eq(n, ns1, ns2, s12, s1, s2), l)

def correlations(i_ang):
	
	ang_low = (pi / 180) * (ang[i_ang] - (ang_res / 2))
	ang_high = (pi / 180) * (ang[i_ang] + (ang_res / 2))

	n = np.zeros(N_jkmpix)	
	c = np.zeros((n_corr, N_jkmpix))
	c_out = np.zeros((n_corr, N_jkmpix))
	#c_gal = 0
	#n_gal = 0
	
	for i in range(N_mpix):
		i_vec = hp.pix2vec(nside, mpix2hpix[i])
		disc_low = hp.query_disc(nside, i_vec, ang_low, inclusive = False)
		disc_high = hp.query_disc(nside, i_vec, ang_high, inclusive = False)
		disc = hpix2mpix[np.setdiff1d(disc_high, disc_low)]
		disc = disc[disc >= 0]
		
		k = hp.vec2pix(nside_jk, i_vec[0], i_vec[1], i_vec[2])
		k = jkhpix2jkmpix[k]

		for l in range(n_corr):		
			dmap1 = dmap[l][0]	
			dmap2 = dmap[l][1]	
			c[l][k] += dmap1[i] * dmap2[disc].sum()

		n[k] += len(disc) 

		#c_gal += dmap[0][0][i] * dmap[0][1][disc].sum()
		#n_gal += len(disc)

	for l in range(n_corr): 
		for k in range(N_jkmpix):
			c_out[l][k] = (c[l].sum() - c[l][k]) / (n.sum() - n[k])

	print "theta = %.2f Done!" % ang[i_ang]
	return c_out

#MAIN...........................................

#Catalog Dictionary structure
#cat['s'] - type of catalog s:spectroscopic p:photometric
#        ['pwd'] - path of the folder where the catalog is located
#        ['f'] - file name of the catalog
#        ['col'] - column index where the values are
#               ['zp'] - array with photometric redshift
#               ['zs'] - array with spectroscopic redshift
#               ['od'] - array with photo-z quality index
#        ['val']['zp', 'zs', 'od', 'm']
#cat['p']
#        ['pwd']
#        ['f']
#        ['col']['zp', 'od', 'ra', 'dec', 'm']
#        ['val']
#               ['zp']
#               ['od']
#               ['ra'] - array with right ascension (in degrees)
#               ['dec'] - array with declination values (in degrees)
#               ['hpix'] - array with healpix index
#               ['mpix'] - array with maskpix index

#Map Dictionary structure
#map - array with equal lenght as id_mask
#   ['stars'] - the number of stars at each mask pix
#   ['gal'] -  the number of galaxies at each mask pix 
#          ['zbinX.XX_X.XX'] - photo-z bin of the map
#                           ['odX.XX'] - od cut of the map
#   ['odds'] - the averaged value of the odds of the galaxies at each mask pix 
#           ['zbinX.XX_X.XX']['odX.XX']

#Angular Correlations Dictionary structure
#corr - array containing the correlation values in each theta bin
#    ['gal-gal'] - maps that are cross-correlated
#               ['zbinsX.XX_X.XX-X.XX_X.XX'] - photo-z bins cross-correlated
#                                           ['odX.XX'] - odds cut applied
#						      ['jkX'] - jacknife region
#    ['odds-odds']['zbinsX.XX_X.XX-X.XX_X.XX']['jkX']
#    ['star-star']['jkX']
#    ['gal-odds']['zbinsX.XX_X.XX-X.XX_X.XX']['odX.XX']['jkX']
#    ['gal-stars']['zbinsX.XX_X.XX-X.XX_X.XX']['odX.XX']['jkX']
#    ['odds-stars']['zbinsX.XX_X.XX-X.XX_X.XX']['jkX']

#hpix2mpix - array of dimension 12 * Nside^2 with -1 in the healpixels that are
#not in the mask and maskpixel id in the pixels of the mask

#mpix2hpix - array of dimension N_mpix that contains the corresponden
#healpix index. 

#hmap = healpix map
#mmap = mask map, the dimension is equal to mpix2hpix + 1. The last element has 
#value 0

#Loading catalogs........................................
for i in cat:
	for j in cat[i]['col']:
		_, cat[i]['val'][j] = np.loadtxt(cat[i]['pwd'] + cat[i]['f'], usecols = (1,cat[i]['col'][j]), unpack = True)

if(cat_list['p'] == True): 
	cat['p']['val']['hpix'] = coord2pix(cat['p']['val']['ra'], cat['p']['val']['dec'], nside)

	#Generating mask map conversion..................
	mask_hmap = hp.read_map(mask_file)
	hpix2mpix = idmaskmap(mask_hmap)
	mpix2hpix = np.nonzero(mask_hmap)[0] 
	N_mpix = len(mpix2hpix)
	print "N_mpix =", N_mpix

	#Generating jacknife map conversion..............
	jk_hmap = hp.ud_grade(mask_hmap, nside_jk)
	jkhpix2jkmpix = idmaskmap(jk_hmap)
	jkmpix2jkhpix = np.nonzero(jk_hmap)[0]
	N_jkmpix = len(jkmpix2jkhpix)
	print "N_jkpix =", len(jk_hmap)
	print "N_jkmpix =", N_jkmpix

	cat['p']['val']['mpix'] = hpix2mpix[cat['p']['val']['hpix']]

for i in cat: 
	#cat[i]['val']['od'] = odrescal(cat[i]['val']['od'], od_range)
	cat[i]['val']['od'] = tls.od_renorm(cat[i]['val']['od'], od_min = od_range[0])

od_cut = tls.od_cut(cat['s']['val']['od'], eff_cut)

if(cat_list['s'] == True): 

	#Generating Nz....................................
	#if not os.path.exists('./Nz'):
	#	os.makedirs('./Nz')
	
	mycolors = tls.spectral(n_od, 0.2,1.)
	for bin in zbin:
		c = 0
		for od in od_cut:
			mask = cut(bin, od, cat['s']['val']['zp'], cat['s']['val']['od'])		

			zs = cat['s']['val']['zs'][mask]
			h, x = np.histogram(zs, bins = 25 , normed = True, range = (0.,1.))
			x = x[:-1] + (x[1] - x[0])/2
		
			np.savetxt(folder_out + 'Nz/Nz' + '_' + zbintag(bin) + '_' + odtag(eff_cut[c]) + '.txt', np.array([x,h]).T, fmt = '%.4f')
		
			H, _ = np.histogram(zs, bins = 20 , range = (.2,.9))
			if(H.sum() > 100): plt.plot(x,h, color = mycolors[c], lw = 2, alpha = 0.8)
			c += 1

	#Ploting Nz......................................
	#if not os.path.exists('./plot'):
	#	os.makedirs('./plot')

	#plt.savefig(folder_out + 'plot/Nz.png')

if(cat_list['p'] == True): 

	#Generating maps.................................
	map = {}
	if((map_list['gal'] == True) or (map_list['odds'] == True)): map['gal'] = {}
	if(map_list['odds'] == True): map['odds'] = {}
	map_generator()
	
	#Computing correlations..........................
	ang = np.arange(ang_min, ang_max, ang_step)
	n_ang = len(ang)

	corr = set_corrdic()
	corr_list = list(paths(corr))
	print corr_list
	n_corr = len(corr_list)

	dmap = np.zeros((n_corr, 2, N_mpix))
	for l in range(n_corr): 
		map1, map2 = get_map(corr_list[l])
		dmap[l][0] = delta(map1) 
		dmap[l][1] = delta(map2) 

	po = Pool() 
	res = po.map(correlations, np.arange(n_ang))
	po.close()
	
	res = np.swapaxes(res, 0,1)
	for l in range(n_corr): write_corr(res[l], corr_list[l])	

	#Corrections..................................
	gal_corr_list = [] 
	for l in range(n_corr): 
		if((corr_list[l][0] == 'gal') and (corr_list[l][1] == 'gal')): gal_corr_list += [corr_list[l]] 
	for s in correction_list:
		for l in gal_corr_list:
			correction(l, s)

	#Writing correlations..........................
	#if not os.path.exists('./corr'):
	#	os.makedirs('./corr')
	
	for l in corr_list: 
		if((l[0] == 'gal') and (l[1] == 'gal')): file_name = folder_out + 'Corr/corr' + nside_tag + corrtag(l) + correctiontag(correction_list) + angtag(ang_min, ang_max) + '.txt'
		else: file_name = folder_out + 'Corr/corr' + nside_tag + corrtag(l) + '.txt'
		val = get_dict(corr,l)	
		mean = np.sum(val, axis = 1) / N_jkmpix 
		cov = np.zeros((n_ang, n_ang))
		for i in range(n_ang):
			for j in range(n_ang):
				for k in range(N_jkmpix):
					cov[i][j] += (val[i][k] - mean[i]) * (val[j][k] - mean[j])
					cov[i][j] = ((float(N_jkmpix) - 1.) / float(N_jkmpix)) * cov[i][j] 
		
		np.savetxt(file_name, np.append([ang,mean], cov, axis =0).T)
