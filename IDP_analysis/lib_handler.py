"""
This script is designed to handle importing. In practice it's not much
better than simply importing all of the packages one-by-one, but I like
the try/except approach with the slightly less grumpy error message



"""


libs = {'numpy':'np','mdtraj':'md','matplotlib':'mpl','matplotlib.pyplot':'plt'}

import importlib, os, sys, shutil, subprocess

for lib in sorted(libs.iteritems()):
	try:
		globals()[lib[0]] = importlib.import_module(lib[0])
		globals()[lib[1]] = globals().pop(lib[0])
	except ImportError:
		print lib, "not found. Did you forget to activate your conda environment?"
		exit()
	# need to set backend before importing plt
	if lib[0] == "matplotlib":
		if 'DISPLAY' not in os.environ:
			mpl.use('Agg')
		font = {'family' : 'normal','weight' : 'normal','size'   : 15}
		mpl.rc('font', **font)
		from matplotlib.colors import LogNorm

