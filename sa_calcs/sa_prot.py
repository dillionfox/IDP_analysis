from utils import np, md, plt, subprocess, os
from prot_geo import prot_geo
from polymer import polymer
from rama import rama
from pca import pca
import kmeans_clusters

class sa_prot(prot_geo,polymer,rama,pca):
	def __init__(self):
		prot_geo.__init__(self)
		polymer.__init__(self)
		rama.__init__(self)
		pca.__init__(self)

		#---Rosetta
		self.scores = []		# list of scores from Rosetta with global index: [ind, score]
		self.ros_frames = -1		# number of "top score" structures to look at. defined with class call (not internal to class)
		self.score_file = ''		# path to Rosetta score file.

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

	def protein_calcs(self,struc,calcs):
		coors = struc.xyz[0] 
		CA_coors = struc.atom_slice(struc.topology.select('name CA'))[0].xyz[0]
		self.nres = struc.n_residues

		if 'Gyr' in calcs:
			L = self.gyration_tensor(coors)
		if 'Rg' in calcs:
			self.compute_Rg(L)
		if 'Asph' in calcs:
			self.compute_Asph(L)
		if 'EED' in calcs:
			self.EED.append(np.linalg.norm(CA_coors[0]-CA_coors[-1]))
		if 'SASA' in calcs:
			SASA = md.shrake_rupley(struc)
			self.SASA.append(SASA.sum(axis=1)[0])
		if 'cmaps' in calcs:
			dist = self.contact_maps(CA_coors)
			self.cmaps.append(dist)
		if 'gcmaps' in calcs:
			self.gremlin_contact_maps(dist)
		if 'SS' in calcs:
			self.SS.append(md.compute_dssp(struc))
		if 'flory' in calcs:
			self.compute_flory(struc,self.nres)
		if 'rama' in calcs:
			self.dihedrals.append(self.compute_phipsi(struc))
		if 'surface_contacts' in calcs:
			self.scmaps.append(self.surface_contacts(struc,SASA))

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

	@staticmethod
	def load_scores(fil): 
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

