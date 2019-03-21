from lib_handler import np, md, plt, subprocess, os
from MASA import MASA
from lipid_analysis import lipid_analysis
from curvature import curvature
from contacts import contacts

class sa_mem(MASA, lipid_analysis, curvature, contacts):
	def __init__(self,outdir,name_mod):
		MASA.__init__(self,outdir,name_mod)
		lipid_analysis.__init__(self,outdir,name_mod)
		curvature.__init__(self,outdir,name_mod)
		contacts.__init__(self,outdir,name_mod)

