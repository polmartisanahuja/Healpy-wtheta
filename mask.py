import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from math import *

#Parameters..................................................
from parameters import *

def mask(hist_pix):
  """Returns an histogram of pixels within the MASK.
  Tricky MASK!!: By now, this is a provisional mask that
  that removes N>5 contiguos pixels with zeros content."""
  
  pix_index = np.nonzero(hist_pix)[0]
  N = len(pix_index)

  #Find zeros_threshold
  print "Consecutive zeros histogram:"
  D = []
  for i in range(1, N): 
    sys.stdout.write("\r\tComplete: %1.1f%%" % (100 * float(i) / float(N)))
    sys.stdout.flush()
    D = np.append(D, pix_index[i] - pix_index[i - 1])
  
  n_bins = 30
  plt.hist(D, bins = n_bins, range = (2, n_bins + 2))
  plt.xlabel("consecutive zeros")
  plt.show()

  zeros_trsh = input("\nzeros_trsh in the MASK:")

  #Re-add zero content pixels due to statistical fluctiations
  print "Re-add statistical zero chains:"
  zeros_index = []
  zeros_pix_index = []
  for i in range(1,N):
    sys.stdout.write("\r\tComplete: %1.1f%%" % (100 * float(i) / float(N)))
    sys.stdout.flush()
    d = pix_index[i] - pix_index[i - 1]
    if ((d < zeros_trsh) & (d > 1)): 
      zeros_index = (np.append(zeros_index, i * np.ones(d - 1))).astype(int)
      zeros_pix_index = (np.append(zeros_pix_index, np.arange(pix_index[i - 1] + 1,  pix_index[i], 1,))).astype(int)
  pix_index = np.insert(pix_index, zeros_index, zeros_pix_index )

  #Create the mask throught pix index
  mask_map = np.zeros(hp.nside2npix(nside))
  for i in pix_index:
    mask_map[i] = 1

  return mask_map 

#MAIN..........................................................

print "******************************************"
print "*                                        *"
print "*           Mask map by P.Marti          *" 
print "*                                        *"
print "******************************************"

print "Reading map..."
hist_pix = hp.read_map(n_map_file)

print "Creating the mask..."
mask_map = mask(hist_pix)

print "\nWriting mask..."
hp.write_map(mask_file, mask_map)
