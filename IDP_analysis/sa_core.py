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

def header():
	"""
	Print header

	"""
	print "                      _________ ______   _______                     "
	print "                      \__   __/(  __  \ (  ____ )                    "
	print "                         ) (   | (  \  )| (    )|                    "
	print "                         | |   | |   ) || (____)|                    "
	print "                         | |   | |   | ||  _____)                    "
	print "                         | |   | |   ) || (                          "
	print "                      ___) (___| (__/  )| )                          "
	print "                      \_______/(______/ |/                           "
	print "                                                                     "
	
	print " _______  _        _______  _              _______ _________ _______  "
	print "(  ___  )( (    /|(  ___  )( \   |\     /|(  ____ \ __   __/(  ____ \ "
	print "| (   ) ||  \  ( || (   ) || (   ( \   / )| (    \/   ) (   | (    \/ "
	print "| (___) ||   \ | || (___) || |    \ (_) / | (_____    | |   | (_____  "
	print "|  ___  || (\ \) ||  ___  || |     \   /  (_____  )   | |   (_____  ) "
	print "| (   ) || | \   || (   ) || |      ) (         ) |   | |         ) | "
	print "| )   ( || )  \  || )   ( || (____/\| |   /\____) |___) (___/\____) | "
	print "|/     \||/    )_)|/     \|(_______/\_/   \_______)\_______/\_______) "
	print
	print
	print "See https://github.com/dillionfox/IDP_analysis for basic information  "
	return None

def usage():
	print "USEAGE: python rosetta_analysis.py ARGS"
	print "ARGS: EITHER a list of pdbs in a file with a .txt extension, or"
	print "a .dcd/.xtc and a .pdb"
