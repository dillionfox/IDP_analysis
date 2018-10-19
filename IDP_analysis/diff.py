import mdtraj as md
import numpy as np
import sys
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def plot_shells(R,D,outdir,name_mod):
	"""
	Very basic plotting function

	"""
	print " ~!~ DIFFUSION ~!~", R, D
	plt.xlabel('nm')
	plt.ylabel('Diffusion Constant')
	plt.plot(R,D)
	plt.axhline(max(D) - 0.1*(max(D)-min(D))) # draw line at 90% mark (arbitrary cutoff)
	plt.subplots_adjust(top=0.93,bottom=0.15,left=0.2,right=0.9,hspace=0.1,wspace=0.1)
	#plt.show()
	plt.savefig(outdir+"diffusion" + name_mod + ".png")
	return 0

def D_shells(fr,fr0,R1,R2,selname="protein"):
	"""
	This function computes the diffusion of water around molecules
	in a discrete shell

	USAGE:
	fr		current frame
	fr0		previous frame
	R1		lower bound
	R2		upper bound
	selname		selection at center of water shells. default is TRG

	RETURN:
	D 		the average distance each water molecule moved (which is related to diffusion)

	"""

	oxygen = fr.topology.select("water and name 0H2")
	sel = fr.topology.select(selname)
	within_R1 = md.compute_neighbors(fr,R1,sel,oxygen)[0]
	within_R2 = md.compute_neighbors(fr,R2,sel,oxygen)[0]
	between_R1_R2 = np.setdiff1d(np.array(within_R2), np.array(within_R1))
	shell_coors = fr.atom_slice(between_R1_R2)[0].xyz[0]
	shell_coors_old = fr0.atom_slice(between_R1_R2)[0].xyz[0]
	norm_list = np.array([np.linalg.norm(i) for i in (shell_coors-shell_coors_old)])
	return np.mean(norm_list[np.where(norm_list < 3)[0]], axis = 0)**2

if __name__ == "__main__":
	dcd = sys.argv[1] ; psf = sys.argv[2]
	traj = md.load_dcd(dcd, top=psf)
	R_list = np.array([0, 0.3, 0.5, 0.6, 1.0])
	diff_data = np.zeros((len(traj)-1,len(R_list),2))
	for i,fr in enumerate(traj):
		if i>0:
			for ri in range(1,len(R_list)):
				diff_data[i-1][ri] = np.array([(R_list[ri]+R_list[ri-1])/2.,\
					np.mean(diff.D_shells(fr,traj[i-1],R_list[ri-1],R_list[ri]),axis=0)])
	#plot_shells(data) # this line needs fixed
