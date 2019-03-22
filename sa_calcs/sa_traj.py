from utils import np, md, plt, os
from diff import diff

class sa_traj(diff):
	def __init__(self,trajname,top):
		#---Structural info
		self.trajname = trajname	# store name of trajectory
		self.top = top			# store topology (pdb, only needed if traj is compressed)
		self.first_frame = 0            # 
		self.last_frame = -1            # NEEDS TO BE FIXED!!
		self.skip_frames = 1		# skip this many frames
		self.nframes = -1		# number of frames in trajectory
		self.nres = -1			# number of residues in structure
		self.protein_analysis = True	# sometimes I run calculations on membrane-only systems
		self.tpr = None			# some membrane analysis calculations (order, density) require gmx make_ndx 
		diff.__init__(self)
		self.traj = None
		self.rmsd = None

	def diffusion(self,fr,calcs):
		"""
		Separate function for diffusion calcs which require protein+water

		"""
		if 'diffusion' in calcs:
			struc = traj[fr_] ; struc_0 = traj[fr_-1] ; N = len(traj)
			# only start on second frame
			if self.diff_data == []:
                        	self.diff_data = np.zeros((N-1,len(self.R_list)))		
			if fr>0:
				for ri in range(1,len(self.R_list)):
					self.diff_data[fr-1][ri] = self.D_shells(struc,struc_0,self.R_list[ri-1],self.R_list[ri])
		return None

	def traj_post_analysis(self,calcs):
		if 'diffusion' in calcs:
			D = np.mean(np.array(self.diff_data).T,axis=1)[1:]
			R = [(self.R_list[i]+self.R_list[i-1])/2. for i in range(1,len(self.R_list))]
			self.plot_shells(R,D,self.outdir,self.name_mod)

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

