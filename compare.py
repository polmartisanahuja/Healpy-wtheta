import numpy as np
import  matplotlib.pyplot as plt
from parameters import *
from enviromental import *

l = 'n_auto_corrected'

for od in odds_cut_list:
	od = "%.2f" % od
	corr_file = "corr/" + l + "_od" + od + ".txt"	
	corr = np.loadtxt(corr_file, unpack = True)
	plt.errorbar(corr[0], corr[1], corr[2],  label = "Odds>" + od)

#plt.ylim(ymax = pow(10,-1), ymin = pow(10,-3))
#plt.xlim(xmax = 3., xmin = 0.)
#plt.semilogy()
#plt.axhline(0)
plt.legend(prop={'size':10}, loc = 'best')
plt.title(l)
plt.xlabel("angle (deg)")
plt.savefig("plot/" + l + ".png")
#plt.show()
plt.close()
