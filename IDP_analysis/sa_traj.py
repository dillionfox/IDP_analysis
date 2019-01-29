import numpy as np
import mdtraj as md
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
mpl.use('Agg')
import matplotlib.pyplot as plt
import os

class sa_traj:
	def __init__(self):
		self.traj = None
		self.rmsd = None

	def traj_calcs(self,traj):
		self.traj = traj
		if 'rmsd' in self.calcs:
			self.rmsd = self.RMSD()
			plt.plot(self.rmsd)
			np.savetxt(self.outdir+"rmsd" + self.name_mod + ".npy",self.rmsd)
			plt.savefig(self.outdir+"rmsd" + self.name_mod + ".png")
		return None

	def RMSD(self,traj):
		"""
		Use MDTraj to compute the RMSD relative to the first frame
	
		"""
		self.rmsd = md.rmsd(traj,traj)
		return None

	@staticmethod
	def read_calibur(fil):
		import shutil
		for line in open(fil):
			line = line.split()
			if len(line) == 0:
				continue
			if line[0] == 'cluster' and line[1] == '=':
				cluster = line[2].split(';')[0]
				center = line[5]
				outdir = '/'.join(center.split('/')[:-1])
				cdir = outdir+'/'+str(cluster)
				if not os.path.exists(cdir):
					os.makedirs(cdir)
				all_ = line[11:]
				for a in all_:
					name = a.split('/')[-1]
					shutil.move(a,cdir+'/'+name)

	@staticmethod
	def calibur(traj,outdir):
		import subprocess
		CALIBUR="/home/dillion/pkg/rosetta/rosetta_bin_linux_2017.08.59291_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/calibur.static.linuxgccrelease"
	
		if not os.path.exists(outdir+"pdbs"):
			os.makedirs(outdir+"pdbs")
		file_list = open(outdir+"pdb_list.txt","w")
		for it,struc in enumerate(traj):
			pdb_name = outdir+"pdbs/"+str(it)+".pdb"
			struc.save(pdb_name)
			file_list.write(pdb_name+"\n")
		file_list.close()
		calibur_out = outdir+"calibur.out"
		FNULL = open(calibur_out,"w")
		subprocess.call([CALIBUR,'-pdb_list',outdir+"pdb_list.txt"], stdout=FNULL, stderr=subprocess.STDOUT)
		read_calibur(calibur_out)

