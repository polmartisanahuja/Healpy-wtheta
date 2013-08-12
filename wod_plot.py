import numpy as np
import matplotlib.pylab as plt
import tools as tls
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from wod_functions import *
from wod_parameters import *

#Plotting correlations........................
plt.figure(1, figsize=(3.5 * n_zbin, 3 * n_zbin))
plt.subplots_adjust(hspace=0,wspace=0)

gs = gridspec.GridSpec(n_zbin, n_zbin)

mycolors = tls.spectral(n_od, 0.2, 1.0)
for i in range(n_zbin):
	for j in range(i, n_zbin):
		ax = plt.subplot(gs[j, i])

		for o in range(n_od):
			corr_file =  file_name = './corr/corr' + nside_tag + '_gal_gal_' + zbintag(zbin[i]) + '_' + zbintag(zbin[j]) + '_' + odtag(od_cut[o]) + correctiontag(correction_list) + angtag(ang_min, ang_max) + '.txt'
			a = np.loadtxt(corr_file, unpack = True)
			theta = a[0]
			corr = a[1]
			cov = a[2:]
			err_corr = np.sqrt(cov.diagonal())

			plt.errorbar(theta, corr, err_corr, color = mycolors[o])
			plt.axhline(y=0.0, color = 'black', linewidth = 1, linestyle = 'dashed')
			if(i == j): plt.title(zbintag(zbin[i]))
			if(i != 0): plt.setp( ax.get_yticklabels(), visible=False)
			else: 
				plt.ylabel('$\omega(\\theta)$')
				plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))	
			if(j != n_zbin - 1): plt.setp( ax.get_xticklabels(), visible=False)
			else: 
				plt.xlabel('$\\theta$(degrees)')
				plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))	
			#plt.yscale('log')
			plt.ylim(ymax = 0.01, ymin = -0.003)
			#plt.ylim(ymax = 0.0065, ymin = 0.)
			#plt.ylim(ymax = 0.105, ymin = -0.002)
			#plt.ylim(ymax = 0.015, ymin = -0.001)
			#plt.xlim(xmax = 7., xmin = 0)
			#plt.ylim(ymax = 0.15, ymin = 0.0001)

plt.savefig('./plot/corr_gal_gal.pdf')	
