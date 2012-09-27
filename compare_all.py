import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from enviromental import *

for l in corr_list:
	title = l 

	for od in odds_cut_list:
		od = "%.2f" % od
		corr_file = "corr/" + l + "_od" + od + ".txt"	
		corr = np.loadtxt(corr_file, unpack = True)
		plt.errorbar(corr[0], corr[1], corr[2],  label = "Odds>" + od)

	plt.axhline(0)
	plt.legend(prop={'size':10}, loc = 'best')
	plt.title(title)
	plt.xlabel("angle (deg)")
	plt.savefig("plot/" + l + ".png")
	plt.close()
