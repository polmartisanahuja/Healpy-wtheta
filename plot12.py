import os
import numpy as np
import matplotlib.pyplot as plt
from parameters12 import *
from enviromental import *

corr_dir_name = './corr/' + dir_name
plot_dir_name = './plot/' + dir_name
if not os.path.exists(plot_dir_name):
	os.makedirs(plot_dir_name)

def plot_n():
	lis = np.array(['n12','n12_corrected'])
	for l in lis:
		i = 0
		for od in Od_cut:
			od_tag = "_od%.2f" % od
			corr_file = corr_dir_name + '/' + cat_name + nside_tag + nside_jk_tag + '_zbin' + z_tag[0] + '_' + z_tag[1] + od_tag + '_' + l +".txt"
			corr = np.loadtxt(corr_file, unpack = True)
			plt.errorbar(corr[0], corr[1], corr[2],  label = "Odds>" + str(od), color = colors[i])

			if(theo_curves):
				theo_corr_file = corr_dir_name + '/' + theo_cat_name + '_zbin' + z_tag[0] + '_' + z_tag[1] + od_tag + '_n12_theo.dat'
				theo_corr = np.loadtxt(theo_corr_file, usecols = (0,3), unpack = True)
				plt.plot(theo_corr[0], theo_corr[1], color = colors[i])
			i += 1
 
		if(y_lim): plt.ylim(ymax = y_lim[1], ymin = y_lim[0])
		if(x_lim): plt.xlim(xmax = x_lim[1], xmin = x_lim[0])
		if(log_scale): plt.semilogy()
		plt.axhline(0, color = 'black', linestyle = 'dashed')
		plt.legend(prop={'size':10}, loc = 'best')
		plt.title(l)
		plt.xlabel("angle (deg)")
		plt.savefig(plot_dir_name + "/" + l + ".png")
		#plt.show()
		plt.close()


def plot_od():
	lis = np.array(['od12','od1','od2'])
	for l in lis:
		title = "Odds correlation" 
		corr_file = corr_dir_name + '/' + cat_name + nside_tag + nside_jk_tag + '_zbin' + z_tag[0] + '_' + z_tag[1] + '_od0.00_' + l +".txt"
		corr = np.loadtxt(corr_file, unpack = True)
		plt.errorbar(corr[0], corr[1], corr[2], label = l)

	plt.axhline(0, color = 'black', linestyle = 'dashed')
	plt.legend(prop={'size':10}, loc = 'best')
	plt.title(title)
	plt.xlabel("angle (deg)")
	plt.savefig(plot_dir_name + "/od.png")
	plt.close()

def plot_nod():
	lis = np.array(['nod1','nod2'])
	for l in lis:
		title = "Gal-odds correlation" 
		for od in Od_cut:
			od_tag = "_od%.2f" % od
			corr_file = corr_dir_name + '/' + cat_name + nside_tag + nside_jk_tag + '_zbin' + z_tag[0] + '_' + z_tag[1] + od_tag + '_' + l +".txt"
			corr = np.loadtxt(corr_file, unpack = True)
			plt.errorbar(corr[0], corr[1], corr[2], label = l + " Odds>" + str(od))

	plt.axhline(0, color = 'black', linestyle = 'dashed')
	plt.legend(prop={'size':10}, loc = 'best')
	plt.xlabel("angle (deg)")
	plt.title(title)
	plt.savefig(plot_dir_name + "/nod.png")
	plt.close()

#Main..........................................
if(pn): plot_n()
if(pod): plot_od()
if(pnod): plot_nod()
