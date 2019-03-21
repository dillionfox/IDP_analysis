from IDP_analysis import sa_core
from IDP_analysis.sa_prot import sa_prot
from IDP_analysis.sa_mem  import sa_mem
from IDP_analysis.sa_traj import sa_traj

"""
			################################
			### Structure Analysis Class ###
			###        Dillion Fox       ###
			###          6/2018          ###
			###        UPenn/ORNL        ###
			################################

This class contains some functions that may be useful for analyzing protein structures.
It can be run several different ways. First, an initial run will compute Rg, SASA, end-to-end distance,
asphericity, average secondary structure per residue, average contact maps, and the flory exponent.

Note that you can pick and choose the tasks you wish to run by specifying the "calcs" variable
upon instantiating the class.

The results can then be analyzed using PCA, which can take in any scalar quantities computed above
and sort them into 'PC' space. The points in PC space are then clustered using k-means. The centroids
of the k-means clusters should be meangingfully distinct, and the variance computed from PCA can tell you 
what distinguishes them.

The structural analyses can then be repeated for each cluster. Another script can be used to organize the
array of plots for each cluster to summarize the structural differences (I have a script that does this on
my github. Look for "rosetta_post.py".

The code can also read in a score.fsc file, which is a direct output file from Rosetta containing estimates
for how "good" each decoy structure is. If this file is read in, the structural analysis can be limited to
the top "N" structures.

This code can also try to fit end-to-end distributions to some fundamental polymer physics models, such
as Gaussian chains or semiflexible polymers.

An example of how I used this script to analyze Rosetta structures can be found at the bottom, under
the "if __name__ == "__main__:" line. Since the trajectories are read by MDTraj, this code is not limited
to dcd's, though some very simple task managing functions would have to be modified to read different 
formats. 

Finally, another script on my github account can be found that uses this class to sequentially analyze 
20 Monte Carlo trajectories of different sequences. See "campari_analyze.py" and "campari_post.py" for
ideas of how this can be used.


"""

# access previously imported libraries
md = sa_core.md
np = sa_core.np
os = sa_core.os

class SA(sa_prot,sa_traj,sa_mem):

	#---Every known calculation
	known_calcs = ['Rg', 'SASA', 'EED', 'Asph', 'rama', 'cmaps', 'PCA', 'gcmaps',\
			'XRg','SS', 'chain', 'score','flory', 'centroids', 'Gyr', \
			'surface_contacts', 'rmsd', 'probe', 'MASA', 'calibur', \
			'diffusion', 'contact_types', 'contact_residues', 'membrane_contacts',\
			'area_per_lipid','persistence_length','membrane_analysis','av_heights',\
			'lipid_phase_transition']

	#---Calculations that depend on each other
	deps = {'Rg':['Gyr'], 'Asph':['Gyr'], 'surface_contacts':['SASA'],\
			'contact_residues':['cmaps'], 'contact_types':['cmaps'],\
			'PCA':['Gyr','Rg','SASA','EED','Asph'],'chain':['Gyr','Rg','EED'],\
			}

	#---Every known post-analysis function
	post_analysis = ['Rg', 'cmaps', 'gcmaps', 'surface_contacts', 'SS' , 'PCA',\
			'flory', 'chain', 'centroids', 'rama', 'MASA', 'diffusion',\
			 'contact_residues', 'contact_types', 'membrane_contacts',\
			'membrane_analysis','av_heights']

	traj_master_list = ['rmsd','calibur','probe','persistence_length','membrane_analysis']

	precalcs_master_list = ['membrane_contacts','membrane_analysis','av_heights','area_per_lipid',\
			'lipid_phase_transition']

	def __init__(self, trajname, top='NULL', name_mod='', outdir='', calcs=[]):
		"""
		Create class objects. Keep most of them empty so they can be populated as the code progresses.

		"""
		#---Inherit classes
		sa_prot.__init__(self)
		sa_traj.__init__(self)
		sa_mem.__init__(self,outdir,name_mod)

		#---Run time options
		self.calcs = np.array(calcs)	# list of calculations to perform at run time
		self.analysis = []		# list of post-processing functions to be run
		self.traj_list = np.array([])	# list of calcs that require whole trajectory
		self.precalcs_list = np.array([])# list of calculations that run before looping through frames
		self.outdir = outdir		# directory to write output to
		self.plot_only = False		# only load data and run analysis
		self.name_mod = name_mod	# several functions are reused, so this name differentiates them
		self.OVERWRITE = 'n'		# overwrite old files? defined with class call
		self.verbose = False		# extra print statements
		self.mode = 'default'		# more than one way to run this code: default, cluster, score (rosetta)

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

		#---Rosetta
		self.scores = []		# list of scores from Rosetta with global index: [ind, score]
		self.ros_frames = -1		# number of "top score" structures to look at. defined with class call (not internal to class)
		self.score_file = ''		# path to Rosetta score file.

	##################################################################
	### Load basic information about the structures being analyzed ###
	##################################################################
	def struc_info(self,struc,nframes):
		"""
		Extract basic information about the structure, i.e. sequence, number of residues

		"""
		try:
			self.get_seq(struc)
		except:
			self.protein_analysis = False
			print "can't read sequence. Did you give it a protein?"
		self.nres = struc.n_residues
		self.nframes = nframes
		return None

	def load_top(self):
		"""
		if PDB list is provided, self.top might need to be supplied. Load first pdb and exit.

		"""
		if self.trajname.split('.')[-1] == 'txt':
			for PDB in open(PDB_list):
				self.top = PDB.split("\n")[0]
				break
		else:
			print "Wrong format. Shouldn't be using this function here..."

		return None

	####################################################
	### Special functions to handle input and output ###
	####################################################
	def check_calcs(self):
		#---If no calculations are specified, run these
		if len(self.calcs) == 0:
			for c in ['Gyr', 'Rg', 'SASA', 'EED', 'Asph', 'PCA']:
				if c not in self.calcs:
					self.calcs = np.append(self.calcs,c)
		for c in self.calcs:
			#---Check for requests that can't be fulfilled
			if c not in SA.known_calcs:
				print c, "is not a known calculation..."
				self.calcs = self.calcs[np.where(self.calcs != c)]
			#---Add dependencies
			if c in SA.deps:
				for ext_dep in SA.deps[c]:
					if ext_dep not in self.calcs:
						self.calcs = np.append(self.calcs,ext_dep)
			#---Add post-processing funcitons
			if c in SA.post_analysis:
				self.analysis.append(c)
			# Pull out functions that don't need to be iterated through individually
			if c in SA.traj_master_list:
				self.traj_list = np.append(self.traj_list,c)
				self.calcs = self.calcs[np.where(self.calcs != c)]
			if c in SA.precalcs_master_list:
				self.precalcs_list = np.append(self.precalcs_list,c)
		analysis_only = ['chain']
		for c in self.calcs:
			if c in analysis_only:
				self.calcs = self.calcs[np.where(self.calcs != c)]
		return None

	def check_input(self):
		"""
		This function determines which calcs need to be run

		"""
		#---Make sure 'outdir' is formatted properly
		try:
			#---Sometimes I accidentally put / on the end twice
			if self.outdir[-2] == '/' and self.outdir[-1] == '/':
				self.outdir = self.outdir[:-1]
			#---Sometimes I over-correct and don't put any
			if self.outdir[-1] != '/':
				self.outdir = self.outdir+'/'
		except:
			if self.outdir == '':
				pass
			else:
				print "something is weird here. outdir =", self.outdir

		#---If outdir doesn't exist, make it
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)

		#---Make sure all calculations are being called properly
		self.check_calcs()

		#---It doesn't make sense to try to overwrite raw data and only plot raw data
		if self.plot_only and self.OVERWRITE:
			print "You can't overwrite raw data and only plot raw data. Continuing, but this will not work."
		#---Diffusion code requires some input
		if 'diffusion' in self.calcs:
			print "reminder: If you're computing the diffusion coefficient from a replica exchange simulation,"
			print "then you must use a continuous trajectory"
			self.R_list = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25,1.5,1.75])	
			# code is not currently set up to compute diffusion from pdbs. Not hard, just not doing it right now.
			if traj_ext == 'txt':
				print "can't compute diffusion constant from pdb's yet. Change the loop structure to fix it"
				exit()
		#---If ros_frames isn't specified, use 100 by default
		if self.mode == 'score' and self.ros_frames == -1:
			self.ros_frames = 100
		return None

	def overwrite(self):
		"""
		If you want to overwrite old data

		"""
		self.OVERWRITE = 'y'
		return None

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

	###########################################################
	# Special run modes: Rosetta Score Mode, Analyze Clusters #
	###########################################################
	def load_scores(self,fil): 
		"""
		Parse Rosetta score file for scores

		"""
		score = [] 
		for line in open(fil): 
			list = line.split() 
			if list[1] == "score": 
				continue 
			score.append(float(list[1])) 
		return np.array(score)

	def analyze_clusters(self):
		"""
		Structures can be characterized by PCA, and the structures can be clustered using k-means in PC-space.
		This function performs a structural analysis on each k-means cluster.
	
		"""
		self.name_mod = '_ros'
		k = np.loadtxt(sa.outdir+'PCA_kmeans_clusters' + sa.name_mod + '.npy')
	
		# Be careful with this. Not all calcs need to be redone (i.e. Rg, SASA, EED, PCA)
		# However, cmaps and SS both need to be recomputed because they were already averaged
		self.calcs = ['flory']
	
		if self.top == "":
			self.load_top()
	
		old_name_mod = self.name_mod
		for cluster in range(int(max(k))+1):
			self.cmaps = []
			self.SS = []
			dcd = md.load_dcd(self.outdir + 'cluster_' + str(cluster) + old_name_mod + '.dcd', self.top)
			self.name_mod = '_cluster_'+str(cluster)
			print "CLUSTER:",cluster
			self.nres = dcd[0].n_residues
			for struc in dcd:
				self.protein_calcs(struc)
			if 'cmaps' in self.calcs:
				self.av_cmaps(self.cmaps,self.nres,self.seq,self.outdir,self.name_mod,"NULL")
			if 'SS' in self.calcs:
				self.av_SS(self.SS)
			if 'flory' in self.calcs:
				self.flory_hist()
		return None

	def run_score_mode(self):
		"""
		This is an alternate way to run the code. This method loads all Rosetta
		scores and picks the structures with the highest score and runs selected
		analyses on them.

		"""
		print "ROSETTA SCORE MODE"
		self.name_mod = '_ros'
		#---Search for score file
		if self.score_file == '':
			success = False
			#---first check the most obvious spot
			if os.path.exists(self.outdir+"score.fsc"):
				scores = self.load_scores(self.outdir+"score.fsc") 
				success = True
			#---maybe it's in a subdirectory
			else:
				for Ind in range(len(self.outdir.split('/'))):
					newpath = ''
					for i in self.outdir.split('/')[:-1*Ind]:
						newpath+=i+'/'
					if os.path.exists(newpath+"score.fsc"): 
						scores = self.load_scores(newpath+"score.fsc") 
						success = True
			#---can't find it using simple methods. Throw an error
			if success == False:
				print "looking for",self.outdir+"score.fsc"
				print "score file could not be located. Exiting..."
				exit()
		#---If name is provided, try to load it. Note: name must be provided using class call (i.e. sa.score_file = '/path/to/score.fsc')
		else:
			print "LOADING SCORE FILE"
			scores = self.load_scores(self.score_file) 
		#---Save all Rosetta ID's with their associated scores
		frame_sel = []
		ind = np.argsort(scores)[:self.ros_frames]
		for i in ind:
			self.scores.append([i, scores[i]])
			frame_sel.append(i)
		#---Save associated files in a separate directory
		if not os.path.exists(self.outdir+"top_score/pdbs"):
			os.makedirs(self.outdir+"top_score/pdbs")
		#---Iterate through pdb's with best scores
		count = 0
		for fr,PDB in enumerate(open(self.trajname)):
			if fr in frame_sel:
				struc = md.load(PDB.split("\n")[0])
				struc.save(self.outdir+"top_score/"+PDB.split("\n")[0].split("/")[-1])
				self.protein_calcs(struc)
				count += 1
		self.top = PDB.split("\n")[0]
		return None

	##################################
	### Run requested calculations ###
	##################################
	def precalcs(self,traj):
		"""
		Things that only need to be computed once at the beginning of a calculation

		'membrane_contacts':self.contacts_precalcs(struc)
		'membrane_analysis':self.make_ndx(self.tpr)
		'av_heights':self.heights_precalcs(struc)

		"""
		print "Running Precalcs", self.precalcs_list
		if 'membrane_contacts' in self.precalcs_list:
			self.contacts_precalcs(traj[0])
		if 'membrane_analysis' in self.precalcs_list:
			self.make_ndx(self.tpr)
		if 'av_heights' in self.precalcs_list:
			self.heights_precalcs(traj[0])
			self.lipid_mesh(traj)
		if 'area_per_lipid' in self.precalcs_list:
			self.apl_precalcs(traj)
		if 'lipid_phase_transition' in self.precalcs_list:
			self.lipid_phase_transition(traj[0])
		return None

	def protein_calcs(self,struc):
		"""
		Run calculations specified in self.calcs. Before running calculation, check
		to make sure it wasn't already done. If it was done before, load the data.

		"""
		coors = struc.xyz[0] 
		CA_coors = struc.atom_slice(struc.topology.select('name CA'))[0].xyz[0]
		self.nres = struc.n_residues

		if 'Gyr' in self.calcs:
			L = self.gyration_tensor(coors)
		if 'Rg' in self.calcs:
			self.compute_Rg(L)
		if 'Asph' in self.calcs:
			self.compute_Asph(L)
		if 'EED' in self.calcs:
			self.EED.append(np.linalg.norm(CA_coors[0]-CA_coors[-1]))
		if 'SASA' in self.calcs:
			SASA = md.shrake_rupley(struc)
			self.SASA.append(SASA.sum(axis=1)[0])
		if 'cmaps' in self.calcs:
			dist = self.contact_maps(CA_coors)
			self.cmaps.append(dist)
		if 'gcmaps' in self.calcs:
			self.gremlin_contact_maps(dist)
		if 'SS' in self.calcs:
			self.SS.append(md.compute_dssp(struc))
		if 'flory' in self.calcs:
			self.compute_flory(struc,self.nres)
		if 'rama' in self.calcs:
			self.dihedrals.append(self.compute_phipsi(struc))
		if 'surface_contacts' in self.calcs:
			self.scmaps.append(self.surface_contacts(struc,SASA))
		return None

	def membrane_calcs(self,struc,fr):
		"""
		Separate function for membrane calcs, which require all atoms, not just protein

		'MASA':self.MASA.append(mem.MASA(struc))
		'membrane_contacts':self.compute_contacts(struc)
		'area_per_lipid':self.compute_area_per_lipid(struc)
		'av_heights':self.heights(struc,fr)

		"""
		if 'MASA' in self.calcs:
			self.compute_MASA(struc,fr)
		if 'membrane_contacts' in self.calcs:
			if not self.plot_only:
				self.compute_contacts(struc)
		if 'area_per_lipid' in self.calcs:
			self.compute_area_per_lipid(struc,fr)
		if 'av_heights' in self.calcs:
			self.heights(struc,fr)
		return None

	def diffusion(self,fr):
		"""
		Separate function for diffusion calcs which require protein+water

		"""
		if 'diffusion' in self.calcs:
			struc = traj[fr_] ; struc_0 = traj[fr_-1] ; N = len(traj)
			# only start on second frame
			if self.diff_data == []:
                        	self.diff_data = np.zeros((N-1,len(self.R_list)))		
			if fr>0:
				for ri in range(1,len(self.R_list)):
					self.diff_data[fr-1][ri] = self.D_shells(struc,struc_0,self.R_list[ri-1],self.R_list[ri])
		return None

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
			self.gmx_density(self.trajname,self.tpr,self.OVERWRITE,self.first_frame,self.last_frame)
			self.gmx_order(self.trajname,self.tpr,self.OVERWRITE,self.first_frame,self.last_frame)
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
		sa_core.header()

		#---Code can currently be run in two modes: default, and 'score' mode
		self.mode = mode

		#---Check to see which calculations need to be run
		self.check_input()
		print self.trajname

		#---Load existing data
		if self.mode == 'default' and self.OVERWRITE == 'n': 
			self.load_data()
		elif self.OVERWRITE == 'y':
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
						self.precalcs(traj)
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
								self.protein_calcs(prot)
							#---Run calculations requiring protein, lipid, and water coordinates
							self.membrane_calcs(struc,fr_)
							#---Special case: protein and water coordinates
							self.diffusion(fr_)
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
			sa_core.usage()
	elif len(sys.argv) == 3:
		if sys.argv[1].split('.')[1] in ['dcd', 'xtc'] and sys.argv[2].split('.')[1] in ['pdb', 'gro']:
			traj = sys.argv[1]
			top = sys.argv[2]
			sa = SA(traj,top,'test','test_traj/',['Rg'])
		else:
			sa_core.usage()
	else:
		sa_core.usage()

	sa.overwrite()
	sa.run()
	#sa.ros_score_sort()
	#kmeans_clusters.write_clusters(sa.k,sa.trajname,sa.top,sa.outdir,sa.name_mod)
	#sa.analyze_clusters()
