import numpy as np
import mdtraj as md

def RMSD(traj):
	"""
	Use MDTraj to compute the RMSD relative to the first frame

	"""
	return md.rmsd(traj)
