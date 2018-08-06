import numpy as np
import mdtraj as md
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
mpl.use('Agg')
import matplotlib.pyplot as plt

def calcs(traj,calcs,name_mod,outdir):
	if 'rmsd' in calcs:
		rmsd = RMSD(traj)
		plt.plot(rmsd)
		np.savetxt(outdir+"rmsd" + name_mod + ".npy",rmsd)
		plt.savefig(outdir+"rmsd" + name_mod + ".png")
	return None

def RMSD(traj):
	"""
	Use MDTraj to compute the RMSD relative to the first frame

	"""
	return md.rmsd(traj,traj)
