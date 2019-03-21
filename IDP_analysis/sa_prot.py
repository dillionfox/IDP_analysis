from lib_handler import np, md, plt, subprocess, os
from prot_geo import prot_geo
from polymer import polymer
from rama import rama
from pca import pca
import kmeans_clusters

class sa_prot(prot_geo,polymer,rama,pca):
	def __init__(self):
		prot_geo.__init__(self)
		polymer.__init__(self)
		rama.__init__(self)
		pca.__init__(self)
