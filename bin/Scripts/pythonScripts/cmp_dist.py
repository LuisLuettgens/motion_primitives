import numpy as np
import sys
dist=float(sys.argv[1])
def f(dist):
	min_dist  = 5107
	max_dist  = 12866
	max_discr = 101
	min_discr = 31
	return int(np.round(min_discr+float(max_discr-min_discr)/(max_dist-min_dist)*(dist-min_dist)))


