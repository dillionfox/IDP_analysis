from sa_calcs import utils
from sa_utils.calc_manager import calc_manager
from sa_calcs.sa_prot import sa_prot
from sa_calcs.sa_mem  import sa_mem
from sa_calcs.sa_traj import sa_traj

"""
			################################
			### Structure Analysis Class ###
			###        Dillion Fox       ###
			###          3/2019          ###
			###        UPenn/ORNL        ###
			################################

This class contains some functions that may be useful for analyzing protein structures and membranes.
You can pick and choose the tasks you wish to run by specifying the "calcs" variable upon instantiating 
the class.

For usage, see the examples in the "example_wrappers" folder.

"""

# access previously imported libraries
md = utils.md
np = utils.np
os = utils.os

class SA(sa_prot,sa_traj,sa_mem,calc_manager):

	def __init__(self, trajname, top='NULL', name_mod='', outdir='', calcs=[]):
		"""
		Create class objects

		"""
		#---Inherit classes
		calc_manager.__init__(self,calcs,outdir,name_mod) 	# make sure the calculations requested make sense and can be done
		sa_traj.__init__(self,trajname,top)			# store trajectory info and some methods that require entire trajectory
		sa_prot.__init__(self)					# main protein attributes and calculations
		sa_mem.__init__(self,outdir,name_mod)			# main membrane attributes and calculations

		#---Rosetta. Will probably write this out soon.
		self.scores = []		# list of scores from Rosetta with global index: [ind, score]
		self.ros_frames = -1		# number of "top score" structures to look at. defined with class call (not internal to class)
		self.score_file = ''		# path to Rosetta score file.

	##################################################################
	### Load basic information about the structures being analyzed ###
	##################################################################
	def load_data(self):
		"""
		Don't compute things twice. Load pre-computed data from previous runs

		"""
		#---Load data
		for c in ['Rg', 'EED', 'Asph', 'SASA', 'cmaps', 'gcmaps', 'rama', 'MASA', 'diffusion', 'flory', 'rmsd', 'membrane_contacts','av_heights']:
			if c in self.calcs and os.path.isfile(self.outdir+c+self.name_mod+file_ext):
				print 'loading data for', c, self.outdir+c+self.name_mod+file_ext
				self.__dict__[c] = np.loadtxt(self.outdir+c+self.name_mod+file_ext)
				self.calcs = self.calcs[np.where(self.calcs != c)]

		#---There's no point in just computing the Gyration Tensor
		if 'Gyr' in self.calcs and ('Rg' not in self.calcs and 'Asph' not in self.calcs):
			self.calcs = self.calcs[np.where(self.calcs != 'Gyr')]

		#---Special cases
		if 'SS' in self.calcs and os.path.isfile(self.outdir+'SS_H'+self.name_mod+file_ext):
			print 'loading data for SS...', self.outdir+'SS_H'+self.name_mod+file_ext
			nres,nframes = np.loadtxt(self.outdir+'SS_H'+self.name_mod+file_ext).shape
			self.SS = np.zeros((3,nres,nframes))
			self.SS[0] = np.loadtxt(self.outdir+'SS_H'+self.name_mod+file_ext)
			self.SS[1] = np.loadtxt(self.outdir+'SS_E'+self.name_mod+file_ext)
			self.SS[2] = np.loadtxt(self.outdir+'SS_C'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'SS')]
		if 'PCA' in self.calcs or 'centroids' in self.calcs:
			# This is a special case. If all of the files already exist, don't bother loading them.
			# The calculation is complete.
			pca1 = self.outdir+'PCA_kmeans_clusters' + self.name_mod + '.npy'
			pca2 = self.outdir+'PCA_kmeans_A' + self.name_mod + '.npy'
			pca3 = self.outdir+'PCA_kmeans_B' + self.name_mod + '.npy'
			if (os.path.isfile(pca1) and os.path.isfile(pca1) and os.path.isfile(pca1)) and 'centroids' not in self.calcs:
				print "skipping PCA calculation..."
				self.calcs = self.calcs[np.where(self.calcs != 'PCA')]
			elif (os.path.isfile(pca1) and os.path.isfile(pca1) and os.path.isfile(pca1)) and 'centroids' in self.calcs:
				print "loading PCA data...", "\n", pca1, "\n", pca2, "\n", pca3
				self.k = np.loadtxt(pca1)
				self.A = np.loadtxt(pca2)
				self.B = np.loadtxt(pca3)
			else:
				pass
		print 'DONE loading data'
		return None

	def write_data(self):
		"""
		This is the core data used by many features of the code. Save it so it doesn't need
		to be recomputed

		"""
		for c in ['Rg', 'EED', 'Asph', 'SASA', 'rama', 'diffusion', 'flory', 'rmsd', 'membrane_contacts','area_per_lipid','av_heights']:
			if c in self.calcs:
				np.savetxt(self.outdir+c+self.name_mod+file_ext,self.__dict__[c])
		return None

	##################################
	### Run requested calculations ###
	##################################
	def traj_calcs(self,traj):
		"""
		Calculations that require all frames at once

		"""
		print "Running", self.traj_list
		if 'rmsd' in self.traj_list:
			self.RMSD(traj)
		if 'calibur' in self.traj_list:
			self.calibur(traj,self.outdir)
		if 'persistence_length' in self.traj_list:
			self.persistence_length(traj)
		if 'membrane_analysis' in self.traj_list:
			self.gmx_density(self.trajname,self.tpr,self.overwrite,self.first_frame,self.last_frame)
			self.gmx_order(self.trajname,self.tpr,self.overwrite,self.first_frame,self.last_frame)
		return None

	def post_process(self):
		"""
		All post-processing functions go here

		"""
		if 'Rg' in self.analysis:
			self.plot_Rg(self.outdir,self.name_mod)
		if 'cmaps' in self.analysis:
			try: av_cmaps = self.av_cmaps(self.cmaps,self.nres,self.seq,self.outdir,self.name_mod,"NULL")
			except: print "CMAPS didnt work"
		if 'gcmaps' in self.analysis:
			try: self.av_cmaps(self.gcmaps,self.nres,self.seq,self.outdir,self.name_mod,"gremlin")
			except: print "grem CMAPS didnt work"
		if 'surface_contacts' in self.analysis:
			self.av_cmaps(self.scmaps,self.nres,self.seq,self.outdir,self.name_mod,"surface")
		if 'SS' in self.analysis:
			self.av_SS(self.SS,self.outdir,self.name_mod) ; return 0
			try: self.av_SS(self.SS,self.outdir,self.name_mod)
			except: print "SS didnt work" ; exit()
		if 'PCA' in self.analysis:
			self.run_PCA(self.EED,self.Rg,self.SASA,self.Asph,self.outdir,self.name_mod,self.mode,self.scores,self.trajname,self.ros_frames,self.analysis)
		if 'flory' in self.analysis:
			self.flory_hist(self.nres,self.outdir,self.name_mod)
		if 'chain' in self.analysis:
			#polymer.persistence_length_2(self.EED,self.outdir,self.name_mod)
			polymer.gaussian_chain(self.EED,self.Rg,self.outdir,self.name_mod)
			polymer.semiflexible_chain(self.EED,self.outdir,self.name_mod)
		if 'centroids' in self.analysis:
			self.cluster_centroids()
		if 'rama' in self.analysis:
			self.rama(self.dihedrals,self.outdir,self.name_mod)
		if 'diffusion' in self.analysis:
			D = np.mean(np.array(self.diff_data).T,axis=1)[1:]
			R = [(self.R_list[i]+self.R_list[i-1])/2. for i in range(1,len(self.R_list))]
			self.plot_shells(R,D,self.outdir,self.name_mod)
		if 'contact_residues' in self.analysis:
			self.contact_residues(av_cmaps,self.seq,self.nframes)
		if 'contact_types' in self.analysis:
			self.contact_types(av_cmaps,self.seq,self.nframes)
		if 'membrane_contacts' in self.analysis:
			self.plot_contact_hist(self.outdir,self.name_mod)
		if 'av_heights' in self.analysis:
			self.normalize_heights()
		if 'membrane_analysis' in self.analysis:
			self.plot_order()
		return None

	#####################
	### Main function ###
	#####################
	def run(self,mode='default'):
	        """
		Runs and handles all function calls. All data is stored in class objects.

	        """
		from timeit import default_timer as timer
		start = timer()

		global file_ext ; file_ext = '_raw.npy'
		global traj_ext ; traj_ext = self.trajname.split('.')[-1]

		#---Print Header
		utils.header()

		#---Code can currently be run in two modes: default, and 'score' mode
		self.mode = mode

		#---Check to see which calculations need to be run
		self.check_input()
		print self.trajname

		#---Load existing data
		if self.mode == 'default' and self.overwrite == 'n': 
			self.load_data()
		elif self.overwrite == 'y':
			print "OVERWRITING OLD DATA!"

		#---Print log of tasks left to complete
		print "calculations left to do:", self.calcs,self.traj_list

		#---Decide if it's necessary to load trajectories/PDBs
		if len(self.calcs) == 0 and len(self.traj_list) == 0: LOAD = False
		else: LOAD = True

		#---Run the code
		if self.mode == 'default' and not self.plot_only:
			#---Set LOAD = False if you just want to post-process existing data
			if LOAD == True:
				print "Loading Trajectory"
				#---Right now the code expects either a list of pdbs (.txt), a .dcd, or a .xtc
				if traj_ext in ['dcd', 'xtc', 'txt']:
					#---XTCs and DCDs will be loaded all at once. Extract basic info about structure
					if traj_ext == 'dcd':
						traj = md.load_dcd(self.trajname, self.top)
						self.struc_info(traj[0],len(traj))
					elif traj_ext == 'xtc':
						traj = md.load_xtc(self.trajname, top=self.top)
						self.struc_info(traj[0],len(traj))
					#---Load names of PDBs to be loaded
					elif traj_ext == 'txt':
						with open(self.trajname) as t:
							self.top = t.readline().rstrip()
							nlines = sum(1 for line in t)
						self.struc_info(md.load(self.top),nlines)
						traj = open(self.trajname)
					#---Pre-calculate some things to avoid computing multiple times
					if traj_ext != 'txt' and len(self.precalcs_list) > 0:
						self.precalcs(traj,self.precalcs_list)
					#---Only load the necessary frames
					if self.last_frame != -1:
						traj = traj[self.first_frame:self.last_frame]
					if self.skip_frames != 1:
						traj = traj[::self.skip_frames]
					#---Calculations done on trajectory all at once
					if len(self.traj_list) > 0:
						self.traj_calcs(traj)
					#---Frame-by-frame calculations
					if len(self.calcs) > 0:
						for fr_,struc in enumerate(traj):
							if self.verbose:
								print "frame", fr_, "of", len(traj)
							#---If .txt, then the structures have to be loaded one-by-one
							if traj_ext == 'txt': struc = md.load(struc.split("\n")[0])
							#---Many calculations only require protein coordinates
							if self.protein_analysis:
								prot = struc.atom_slice(struc.topology.select('protein'))[0]
								#---Run calculations that only require protein coordinates
								self.protein_calcs(prot,self.calcs)
							#---Run calculations requiring protein, lipid, and water coordinates
							self.membrane_calcs(struc,fr_,self.calcs)
							#---Special case: protein and water coordinates
							self.diffusion(fr_,self.calcs)
				#---Write data
				self.write_data()

		#---Code can be run in special mode where it references a Rosetta score file and only computes statistics on top N structures
		elif self.mode == 'score':
			self.run_score_mode() # I copy and pasted the code above without modification. If it doesn't work, move it back.

		#---Run post-processing functions, i.e. plotting, etc.
		self.post_process()

		print "Total execution time:", timer()-start
		return None

###############################
### Obsolete. Needs updated ###
###############################
if __name__ == "__main__":
	import sys
	if len(sys.argv) == 2:
		if sys.argv[1].split('.')[1] == 'txt':
			PDB_list = sys.argv[1]
			sa = SA(PDB_list,'','_test','test',['MASA'])
		else:
			utils.usage()
	elif len(sys.argv) == 3:
		if sys.argv[1].split('.')[1] in ['dcd', 'xtc'] and sys.argv[2].split('.')[1] in ['pdb', 'gro']:
			traj = sys.argv[1]
			top = sys.argv[2]
			sa = SA(traj,top,'test','test_traj/',['Rg'])
		else:
			utils.usage()
	else:
		utils.usage()

	sa.overwrite()
	sa.run()
	#sa.ros_score_sort()
	#kmeans_clusters.write_clusters(sa.k,sa.trajname,sa.top,sa.outdir,sa.name_mod)
	#sa.analyze_clusters()
