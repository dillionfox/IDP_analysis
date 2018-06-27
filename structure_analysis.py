import numpy as np
import mdtraj as md
import os
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

################################
### Structure Analysis Class ###
###        Dillion Fox       ###
###          6/2018          ###
###        UPenn/ORNL        ###
################################

class SA:
	"""
	This class contains some functions that may be useful for analyzing protein structures
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

	def __init__(self, trajname, top='NULL', name_mod='', outdir='', calcs=[]):
		"""
		Create class objects. Keep most of them empty so they can be filled as the code progresses.

		"""
		self.trajname = trajname	# store name of trajectory
		self.top = top			# store topology (pdb, only needed if traj is compressed)
		self.name_mod = name_mod	# several functions are reused, so this name differentiates them
		self.EED = [] 			# end-to-end distance
		self.Rg = [] 			# radius of gyration
		self.XRg = [] 			# predicted x-ray radius of gyration (from crysol)
		self.SASA = [] 			# solvent accessible surface area
		self.Asph = [] 			# asphericity 
		self.SS = []			# secondary structure (per residue)
		self.cmaps = []			# contact maps
		self.gcmaps = []		# contact maps for specific residues
		self.scmaps = []		# contact maps for residues on surface only
		self.dihedrals = []		# store dihedrals to generate Ramachandran plot
		self.calcs = calcs		# list of calculations to perform at run time
		self.outdir = outdir		# directory to write output to
		self.fex = []			# Flory Exponent
		self.k = []			# k-means clusters labels
		self.B = []			# matrix in PC space
		self.scores = []		# list of scores from Rosetta with global index: [ind, score]
		self.ros_calcs = []		# list of calcs for rosetta score analysis
		self.mode = 'default'		# more than one way to run this code: default, cluster, score (rosetta)
		self.resnames = []		# this is populated in the 'surface_conacts' procedure. list of residue names for sequence
		self.ind_list = []		# rearranged input axes for PCA
		self.ros_frames = -1		# number of "top score" structures to look at. defined with class call (not internal to class)
		self.OVERWRITE = 'n'		# overwrite old files? defined with class call
		self.fex2 = []			# Flory Exponent
		self.score_file = ''		# path to Rosetta score file.

	#####################################
	#####      Core Functions      ######
	#####################################
	def gyration_tensor(self,coors):
		"""
		3x3 tensor. It's just an outerproduct of r and r. COORDINATES MUST BE CENETERED/

		S_nm = 1/N Sum_i { r(i)_n r(i)_m }

		"""
                coors -= np.mean(coors,axis=0)
                S = np.sum([np.einsum('i,j->ij',c,c) for c in coors], axis=0)/(len(coors))
                L = np.linalg.eig(S)[0]
		return np.sqrt(L)

	def compute_Rg(self,L):
		"""
		L's are sqrt(eigenvalues of gyration tensor)
		Rg^2 = L1^2 + L2^2 + L3^2

		"""
		self.Rg.append(np.sqrt(np.dot(L,L)))
                return None

	def compute_Asph(self,L):
	        """
		Compute the Asphericitiy
		Eq From: Simulation Analysis of the Temperature Dependence of Lignin Structure and Dynamics

	        """
		delta = ((L[0]-L[1])**2+(L[1]-L[2])**2+(L[0]-L[2])**2)/(2*sum(L)**2)
		self.Asph.append(delta)
		return None

	def av_SS(self):
		"""
		Average over all Seconday Structure

		"""

		nframes = len(self.SS) ; nres = len(self.SS[0][0])
		H = np.zeros((nres,nframes)) # alpha helix
		E = np.zeros((nres,nframes)) # extended beta sheet
		C = np.zeros((nres,nframes)) # unstructured coil

		fr = 0
		for ss in self.SS:
			h = np.where(ss[0] == 'H')[0]
			for hh in h:
				H[hh][fr] = 1

			e = np.where(ss[0] == 'E')[0]
			for ee in e:
				E[ee][fr] = 1

			c = np.where(ss[0] == 'C')[0]
			for cc in c:
				C[cc][fr] = 1
			fr+=1

		np.savetxt(self.outdir+"SS_H" + self.name_mod + "_raw.npy",H)
		np.savetxt(self.outdir+"SS_E" + self.name_mod + "_raw.npy",E)
		np.savetxt(self.outdir+"SS_C" + self.name_mod + "_raw.npy",C)

		self.plot_SS(H,E,C,nframes,nres)
		return None

	def plot_SS(self,H,E,C,nframes,nres):
		"""
		Plot average secondary structure per residue with error bars

		"""
		H_av = np.zeros(nres) ; H_err = np.zeros(nres) 
		E_av = np.zeros(nres) ; E_err = np.zeros(nres) 
		C_av = np.zeros(nres) ; C_err = np.zeros(nres) 
		for ind in range(nres):
			H_av[ind] = sum(H[ind,:])/nframes ; H_err[ind] = np.std(H[ind,:])
			E_av[ind] = sum(E[ind,:])/nframes ; E_err[ind] = np.std(E[ind,:])
			C_av[ind] = sum(C[ind,:])/nframes ; C_err[ind] = np.std(C[ind,:])

		plt.clf()
		plt.plot(H_av,label='Helix',c='r') ;      plt.errorbar(range(nres), H_av,yerr=H_err/(np.sqrt(nframes)-1), fmt='o',color='r')
		plt.plot(E_av,label='Beta Sheet',c='b') ; plt.errorbar(range(nres), E_av,yerr=E_err/(np.sqrt(nframes)-1), fmt='o',color='b')
		plt.plot(C_av,label='Coil',c='k') ;       plt.errorbar(range(nres), C_av,yerr=C_err/(np.sqrt(nframes)-1), fmt='o',color='k')
		plt.legend(loc=1)
		plt.savefig(self.outdir+"SS" + self.name_mod + ".png")
		np.savetxt(self.outdir+"SS_H_av" + self.name_mod + ".npy",[H_av,H_err])
		np.savetxt(self.outdir+"SS_E_av" + self.name_mod + ".npy",[E_av,E_err])
		np.savetxt(self.outdir+"SS_C_av" + self.name_mod + ".npy",[C_av,C_err])
		return None

	def contact_maps(self,coors):
		"""
		MDAnalysis takes all of the fun out of this one, but the mdad routines are all done within a KD Tree so
		it's faster than anything I would write.

		"""
		import MDAnalysis.analysis.distances as mdad

		dist = mdad.distance_array(coors, coors)
		self.cmaps.append(dist)
		return dist

	def gremlin_contact_maps(self,dist):
		"""
		Specifically search for contacts that were used as constraints derived by GREMLIN

		"""

		print dist
		contact_cutoff = 10
		gremlin = [[6,13],[9,22],[15,19],[14,18],[3,11],[34,40],[3,23],[36,40],[9,13],[25,28],[12,15],[11,23], \
			[26,35],[12,18],[2,5],[17,21],[14,22],[6,9],[41,44],[15,18],[25,30],[9,16],[29,32],[30,33],[6,16]] 

		contacts = np.zeros(dist.shape)
		for n in range(dist.shape[0]):
			for m in range(dist.shape[1]):
				if dist[n][m] < contact_cutoff and (gremlin.count([n,m]) == 1 or gremlin.count([m,n]) == 1):
					contacts[n][m] = 1

		self.gcmaps.append(contacts)
		return None

	def surface_contacts(self,struc,SASA):
		"""
		SASA per residue from SASA per atom

		"""
		plot_single = False

		cmap_cutoff = 5.0 # nm
		sasa_cutoff = 0.5
		# identify surface residues
		S = []
		N = np.array(range(struc.n_residues))
		for r in N:
			ind = struc.topology.select('resid '+str(r))
			SASA_per_res = np.sum(SASA[0][ind[0]:ind[-1]],axis=0)
			if SASA_per_res > sasa_cutoff:
				S.append(SASA_per_res)
			else:
				S.append(0.0)

		#import ipdb; ipdb.set_trace()
		surf_res = N[np.where(np.array(S) > 0)]
		cmap = self.contact_maps(struc.atom_slice(struc.topology.select('name CA'))[0].xyz[0])

		for n in N:
			for m in N:
				if (n not in surf_res or m not in surf_res) or n < m:
					cmap[n][m] = 0
				elif cmap[n][m] < cmap_cutoff:
					cmap[n][m] = 1

		if plot_single == True:
			hydrophobic = ['GLY', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'CYS']
			hydrophilic = ['SER', 'THR', 'ASN', 'GLN', 'HIS']
			poscharge =   ['ARG', 'LYS']
			negcharge =   ['ASP', 'GLU']

			for it,rn in enumerate(resnames):
				if rn in hydrophobic:
					plt.axhline(y=it,c='grey',linewidth=2)
					plt.axvline(x=it,c='grey',linewidth=2)
				elif rn in hydrophilic:
					plt.axhline(y=it,c='g',linewidth=     2)
					plt.axvline(x=it,c='g',linewidth=     2)
				elif rn in poscharge:
					plt.axhline(y=it,c='b',linewidth=     2)
					plt.axvline(x=it,c='b',linewidth=     2)
				elif rn in negcharge:
					plt.axhline(y=it,c='r',linewidth=     2)
					plt.axvline(x=it,c='r',linewidth=     2)
				else:
					print "unknown restype:", rn
			scat = np.array(zip(np.where(np.array(cmap) > 0)[0], np.where(np.array(cmap) > 0)[1])).T
			plt.scatter(scat[0], scat[1],marker='s',c='k',s=90)
			#plt.imshow(cmap,cmap='Greys')
			plt.show()

		self.scmaps.append(cmap)

		return None

	def av_cmaps(self,cmaps,mtype="NULL"):
		"""
		Average contact maps. This function has the contact cutoff hard-coded. This
		function works for GREMLIN contact maps too.

		"""
		plt.clf()
		nframes = len(cmaps) 
		av = np.zeros((cmaps[0].shape))

		# save cmaps to npy file. Data must first be reshaped.
		if mtype == "NULL":
			cmaps = np.array(cmaps)
			resh = cmaps.reshape(cmaps.shape[0],cmaps.shape[1]*cmaps.shape[2])
			np.savetxt(self.outdir+"CMAPS" + self.name_mod + "_raw.npy",resh)

			for i in range(cmaps[0].shape[0]):
				for j in range(cmaps[0].shape[1]):
					# for each element of the matrix
					if j > i: # don't compute things twice
						l = []
						for c in cmaps:
							# for each map, determine if there was a contact at that position
							if c[i][j] < 0.7: # nm
								l.append(1)
							else:
								l.append(0)
						av[i][j] = np.std(l)/(np.sqrt(nframes)-1)
						av[j][i] = np.mean(l)
					# dont consider contacts from neighbors
					if i == j or abs(i-j) <= 2:
						av[i][j] = 0
						av[j][i] = 0

		else:
			for m in range(self.nres):
				for n in range(self.nres):
					for fr in range(nframes):
						av[n][m] += self.cmaps[fr][m][n]

			av/=nframes

		fig, ax = plt.subplots()

		plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
		if mtype == "gremlin":
			im = ax.imshow(av, cmap='PuBu')
			cbar = fig.colorbar(im)
			ax.set_title("Average Contact Maps from Rosetta+Gremlin Output")
			plt.savefig(self.outdir+"gremlin_compare_CMAPS" + self.name_mod + ".png")
			np.savetxt(self.outdir+"gremlin_CMAPS" + self.name_mod + "_av.npy",av)
		elif mtype == "surface":
			hydrophobic = ['GLY', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'CYS']
			hydrophilic = ['SER', 'THR', 'ASN', 'GLN', 'HIS']
			poscharge =   ['ARG', 'LYS']
			negcharge =   ['ASP', 'GLU']

			for it,rn in enumerate(self.resnames):
				if rn in hydrophobic:
					plt.axhline(y=it,c='yellow',linewidth=1.5)
					plt.axvline(x=it,c='yellow',linewidth=1.5)
				elif rn in hydrophilic:
					plt.axhline(y=it,c='g',linewidth=1.5)
					plt.axvline(x=it,c='g',linewidth=1.5)
				elif rn in poscharge:
					plt.axhline(y=it,c='b',linewidth=1.5)
					plt.axvline(x=it,c='b',linewidth=1.5)
				elif rn in negcharge:
					plt.axhline(y=it,c='r',linewidth=1.5)
					plt.axvline(x=it,c='r',linewidth=1.5)
				else:
					print "unknown restype:", rn
			ax.set_title("Average Contact Maps of Surface Residues")
			im = ax.imshow(av, cmap='Greys')
			cbar = fig.colorbar(im)
			plt.savefig(self.outdir+"surface_CMAPS" + self.name_mod + ".png")
			np.savetxt(self.outdir+"surface_CMAPS" + self.name_mod + "_av.npy",av)
		else:
			im = ax.imshow(av)
			cbar = fig.colorbar(im)
			ax.set_title("Average Contact Maps")
			plt.savefig(self.outdir+"CMAPS" + self.name_mod + ".png")
			np.savetxt(self.outdir+"CMAPS" + self.name_mod + "_av.npy",av)
			
		return None

	def compute_XRg(self, PDB):
		"""
		X-Ray experiments return higher values of Rg because they include some of the water in the shell. The EMBL
		Program "Crysol" computes a theoretical scattering curve for a protein and returns the Rg.

		"""
		import subprocess

		f = PDB.split('.')[0]
		FNULL = open(os.devnull, 'w')
		subprocess.call(['crysol',f+'.pdb'], stdout=FNULL, stderr=subprocess.STDOUT)
		for line in open(f+'00.log'):
			if "Rg ( Atoms - Excluded volume + Shell ) ................. :" in line:
				self.XRg.append(float(line.split(' : ')[1]))
				os.remove(f+'00.log') ; os.remove(f+'00.alm') ; os.remove(f+'00.int')
		return None

	def scatterplot(self,X,Y,xlabel,ylabel,filename):
		"""
		General scatter plot function

		"""
		plt.clf()
		plt.scatter(X,Y)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.savefig(self.outdir+filename + self.name_mod + ".png")
		np.savetxt(self.outdir+ xlabel + self.name_mod + ".npy",x)
		np.savetxt(self.outdir+ ylabel + self.name_mod + ".npy",y)
		return None

	#####################################
	######           PCA           ######
	#####################################
	def plot_Kmeans_PCA(self, A, B, ind_list, PC_labels, kmeansoutput):
		"""
		Make plot in PC space with points colored by k-means

		"""
		plt.clf()
		plt.figure('Top 2 PCAs colored by K-Means Clustering')
		compare_2 = 'n'
		#---Plot in PCA Space
		if compare_2 == 'y':
			plt.scatter(B[:,0][:50], B[:,1][:50], c=kmeansoutput.labels_[:50],edgecolor='k')
			plt.scatter(B[:,0][50:], B[:,1][50:], c=kmeansoutput.labels_[50:],marker="s")
		else:
			plt.scatter(B[:,0], B[:,1], c=kmeansoutput.labels_)

		plt.xlabel('PC-1')
		plt.ylabel('PC-2')
		plt.xlim([-1,1])
		plt.ylim([-1,0.98])
		plt.title('Top 2 PCAs colored by K-Means Clustering')
		plt.subplots_adjust(top=0.9,bottom=0.15,left=0.18,right=0.85,hspace=0.2,wspace=0.2)
		plt.savefig(self.outdir+'PCA_kmeans' + self.name_mod + '.png')
		if self.name_mod == '_ros_score':
			ind = np.array(self.scores).T[0]
			global_labels = [int(ind[l]) for l in range(self.ros_frames)]
			#pdbs = [PDB.split("\n")[0] for PDB in open(self.trajname)]
			#label_names = [pdbs[l] for l in global_labels]
			np.savetxt(self.outdir+'cluster_labels'+self.name_mod+ '.txt',global_labels,fmt="%s")
		np.savetxt(self.outdir+'PCA_kmeans_clusters' + self.name_mod + '.npy',np.array(kmeansoutput.labels_))
		np.savetxt(self.outdir+'PCA_kmeans_A' + self.name_mod + '.npy',np.array(A))
		np.savetxt(self.outdir+'PCA_kmeans_B' + self.name_mod + '.npy',np.array(B))

	def extract_PCA(self, pca, PC_labels):
		"""
		Figure out which components are the most important

		"""
		elbow_cutoff = 0.85
		if self.mode == 'score':
			elbow_cutoff = 0.8

		print "PCA Variance:", pca.explained_variance_ratio_.cumsum()
		i = 1
		for pc in pca.explained_variance_ratio_.cumsum():
			if float(pc) > elbow_cutoff:
				break
			i+=1
		#print "By the Elbow Rule, we will use", i, "pc's"
		nPCA = i

		comps = pca.components_[0]
		print comps
		ind_list = []
		for c in range(nPCA):
			ind = np.where(comps == max(comps))[0][0]
			ind_list.append(ind)
			comps[ind] = -1

		fil = open(self.outdir+'PCA'+self.name_mod+ '.txt','w')
		fil.write("Important Components:")
		for label in range(len(ind_list)):
			fil.write(PC_labels[ind_list[label]]+",")

		return ind_list

	def compute_Kmeans(self,B):
		"""
		Cluster the points in PC space

		"""
		from sklearn.cluster import KMeans
		elbow_cutoff = 0.7
		if self.mode == 'score':
			elbow_cutoff = 0.75
		Nc = range(1, 20)
		kmeans = [KMeans(n_clusters=i) for i in Nc]
		score = [kmeans[i].fit(B).score(B) for i in range(len(kmeans))]

		min_s = min(score)
		max_s = max(score)
		norm_score = [(s-min_s)/(max_s-min_s) for s in score]
		j = 1

		for s in norm_score:
			if s > elbow_cutoff:
				break
			j+=1

		#print "By the Elbow Rule, we will use", j, "Clusters for K-Means"

		#---Plot Elbow Curve
		#plt.clf()
		#plt.plot(Nc,score)
		#plt.xlabel('Number of Clusters')
		#plt.ylabel('Score')
		#plt.title('Elbow Curve')
		#plt.savefig('kmeans-elbow_curve.png')

		return j

	def norm_shift(self, vec):
		"""
		Force all data to extend from -1 to 1

		"""
		vec = np.array(vec)		# [a, b]
		vec -= min(vec)			# [0, b-a]
		vec /= (max(vec) - min(vec)) 	# [0, 1]
		vec *= 2			# [0, 2]
		vec -= 1			# [-1, 1]
		return vec

	def compute_PCA(self, A):
	        """
		perform Principle Component Analysis
		Borrowed from: https://machinelearningmastery.com/calculate-principal-component-analysis-scratch-python/
		NOT CURRENTLY BEING USED. 

	        """
		M = np.mean(A.T, axis=1)
		C = A - M
		V = np.cov(C.T)
		vectors = np.linalg.eig(V)[1]
		P = vectors.T.dot(np.transpose(C))

		# to run this function, move the following lines to the "run" function
		#---My primitive implementation of PCA
		#M_PCA = self.compute_PCA(A)
		#print "mine:\n",  np.real(M_PCA).T
		#print "sklearn\n", B
		return P

	def run_PCA(self):
		"""
		Main function for running PCA. This function also uses k-means clustering on the PC's.

		"""
		from sklearn.decomposition import PCA
		from sklearn.cluster import KMeans

		# normalize and shift all vectors to be centered around zero
		norm_EED =self.norm_shift(self.EED)
		norm_Rg = self.norm_shift(self.Rg)
		norm_SASA = self.norm_shift(self.SASA) 
		norm_Asph = self.norm_shift(self.Asph)

		#---Prepare array containing all SCALED data
		A = np.array([np.array(norm_EED), np.array(norm_Rg), np.array(norm_SASA), np.array(norm_Asph)]).T

		#---Prepare array containing all UNSCALED data (for plotting purposes later)
		Au = np.array([np.array(self.EED), np.array(self.Rg), np.array(self.SASA), np.array(self.Asph)]).T

		PC_labels = ['End-to-End Distance', 'Radius of Gyration', 'SASA', 'Asphericity']

		#---Do PCA 
		pca = PCA(len(PC_labels))
		pca.fit(A)
		self.B = pca.transform(A)

		#---ind_list contains the important components
		ind_list = self.extract_PCA(pca,PC_labels)
		# FIXFIXFIX
		self.ind_list = ind_list

		#---Do initial K-means clustering to determine number of clusters
		nK = self.compute_Kmeans(self.B)

		#---Use optimum number of clusters for k-means
		kmeans=KMeans(n_clusters=nK)
		if len(ind_list) > 1:
			kmeansoutput=kmeans.fit(np.array([self.B[:,0],self.B[:,1]]).T)

			#---Plot top 2 PCA clusters colored by kmeans
			self.plot_Kmeans_PCA(Au,self.B,ind_list,PC_labels,kmeansoutput)
			self.k = np.array(kmeansoutput.labels_)

			self.cluster_centroids() 

			return self.k
		else:
			#---Don't bother continuing if there's only one PC. Not much to do in this case...
			print "Only one PC:", PC_labels[ind_list[0]]
			print "Not generating plot, but will generate dummy files"
			np.savetxt(self.outdir+'PCA_kmeans_clusters' + self.name_mod + '.npy',np.array([]))
			np.savetxt(self.outdir+'PCA_kmeans_A' + self.name_mod + '.npy',np.array(Au))
			np.savetxt(self.outdir+'PCA_kmeans_B' + self.name_mod + '.npy',np.array(self.B))

		return None

	#####################################
	##### K-Means Cluster Analysis ######
	#####################################
	def write_clusters(self):
		"""
		This function writes a new dcd trajectory for each k-means cluster defined by PCA functions.
		This might be useful if you want to perform structural analyses on individual clusters.

		"""
		print 'loaded', sa.outdir+'PCA_kmeans_clusters' + sa.name_mod + '.npy'
		k = np.loadtxt(sa.outdir+'PCA_kmeans_clusters' + sa.name_mod + '.npy')

		# Check to see if cluster trajectories already exist
		count = 0
		for cluster in range(int(max(k))+1):
			if os.path.isfile(self.outdir + 'cluster_' + str(cluster) + self.name_mod + '.dcd'):
				count+=1

		if count == int(max(k)+1):
			print "cluster trajectories already exist"
			return None

		#---DCD
		if self.trajname.split('.')[-1] == 'dcd':
			dcd = md.load_dcd(self.trajname, self.top)
			for cluster in range(int(max(k))+1):
				with md.formats.DCDTrajectoryFile(self.outdir + 'cluster_' + str(cluster) + self.name_mod + '.dcd', 'w') as f: 
					for frames in np.where(k==cluster):
						for fr in frames:
							f.write(dcd[fr].xyz*10)

		#---PDB
		elif self.trajname.split('.')[-1] == 'txt':
			struc = []
			for PDB in open(PDB_list):
				struc.append(md.load(PDB.split("\n")[0]))

			for cluster in range(int(max(k))+1):
				with md.formats.DCDTrajectoryFile(self.outdir + 'cluster_' + str(cluster) + self.name_mod + '.dcd', 'w') as f: 
					for frames in np.where(k==cluster):
						for fr in frames:
							f.write(struc[fr].xyz*10)
					print "done with cluster", cluster
		return None

	def analyze_clusters(self):
		"""
		This is the "run" function for the cluster analysis. It calls and manages the core functions.

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
				self.calculations(struc)

			if 'cmaps' in self.calcs:
				self.av_cmaps(self.cmaps)
			if 'SS' in self.calcs:
				self.av_SS()
			if 'flory' in self.calcs:
				self.fex_hist()
		return None

	def extract_inds(self,fil):
		"""
		If PCA was previously computed but is being used now (i.e. for cluster_centroids), the
		labels of the data with the largest variance can be re-loaded

		"""
		PC_labels = {'End-to-End Distance':0, 'X-Ray Radius of Gyration':1, 'Radius of Gyration':1, 'SASA':2, 'Asphericity':3, '':4}
		for line in open(fil,"r"):
			comps = line.split(":")[1]
			comps = comps.split(",")
			if len(comps) > 1:
				comp1 = comps[0]	
				comp2 = comps[1]
				return [PC_labels[comp1], PC_labels[comp2]]
			else:
				return [comp1]

	def cluster_centroids(self):
		"""
		Identify the centroids of each k-means cluster in PC space

		"""

		replot_PCA = 'y'

		def find_nearest(array, centroid):
			array = np.array(array) ; centroid = np.array(centroid)
			min = np.linalg.norm(array[0] - centroid)+1
			ind = 0
			for el in array:
				dist = np.linalg.norm(el-centroid)
				if dist < min:
					min = dist
					closest = el
					c_i = ind
				ind+=1
			return [c_i, closest]
		
		if 'PCA' not in self.calcs:
			ind_list = self.extract_inds(self.outdir+'PCA'+self.name_mod+ '.txt')

		N = int(max(self.k))
		centroids = []
		labels = []
		coors = []
		B = np.array([self.B[:,0], self.B[:,1]]).T
		for n in range(N+1):
			b = B[np.where(self.k == n)]
			val = np.mean(b,axis=0)
			centroids.append(val)
			[b_i, centroid] = find_nearest(b,val)
			B_i = np.where(B == centroid)[0][0]
			labels.append(B_i)
			coors.append(b[b_i])
		if self.mode == 'score':
			ind = np.array(self.scores).T[0]
			for l in labels:
				print "len ind", len(ind), "l", l
			global_labels = [int(ind[l]) for l in labels]
			pdbs = [PDB.split("\n")[0] for PDB in open(self.trajname)]
			label_names = [pdbs[l] for l in global_labels]
			np.savetxt(self.outdir+'cluster_centroid_labels'+self.name_mod+ '.txt',label_names,fmt="%s")
		else:
			np.savetxt(self.outdir+'cluster_centroid_labels'+self.name_mod+ '.txt',labels, fmt='%i')

		if replot_PCA in ['y', 'yes', 'on']:
			centroids = np.array(centroids).T
			coors = np.array(coors).T
			plt.scatter(self.B[:,0], self.B[:,1], c=self.k,s=30)
			plt.scatter(centroids[0], centroids[1], c='k',s=30)
			plt.scatter(coors[0], coors[1], c='r',s=30)
			for i in range(N+1):
				plt.plot([centroids[0][i], coors[0][i]], [centroids[1][i], coors[1][i]])
			print "NUMBER", N
			plt.xlim([-2,2])
			plt.ylim([-2,2])
			plt.xlabel('PC-1')
			plt.ylabel('PC-2')
			#plt.show()
			plt.savefig(self.outdir+'cluster_centroid_labels'+self.name_mod+ '.png')
		return None

	#####################################
	#######  Load Rosetta Scores  #######
	#####################################
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

	#####################################
	######### Ramachandran Plot ########
	#####################################
	def compute_phipsi(self,struc,one=False):
		"""
		compute dihedrals for Ramachandran plot

		"""
		phi = md.compute_phi(struc)[1][0]
		psi = md.compute_psi(struc)[1][0]
		dihedrals = np.array([phi,psi]).T
		if one == False:
			self.dihedrals.append(dihedrals)
			return None
		else:
			return dihedrals

	def rama(self):
		"""
		Compute Ramachandran angles

		"""
		self.dihedrals = np.array(self.dihedrals)
		a,b,c = self.dihedrals.shape
		angles = self.dihedrals.reshape(a*b,c).T

		fig, ax = plt.subplots()
		plot = plt.hist2d(angles[0], angles[1], bins = 30, norm=LogNorm(), cmap='hot')
		ax.set_xlabel("phi")
		ax.set_ylabel("psi")
		cbar = plt.colorbar()
		cbar.set_label("number of occurances")
		#plt.show()
		plt.savefig(self.outdir+"RAMA" + self.name_mod + ".png")
		np.savetxt(self.outdir+"RAMA_all" + self.name_mod + ".npy",angles)
		return None

	def plot_struc_rama(self,struc,name):
		"""
		Make Ramachandran heatmap

		"""
		di = np.array(self.compute_phipsi(struc,True))
		self.dihedrals = np.array(self.dihedrals)

		fig, ax = plt.subplots()
		plot = plt.hist2d(self.dihedrals[0], self.dihedrals[1], bins = 30, norm=LogNorm(), cmap='hot')
		ax.set_xlabel("phi")
		ax.set_ylabel("psi")
		cbar = plt.colorbar()
		cbar.set_label("number of occurances")
		for i, txt in enumerate(range(1,len(di.T[0]))):
			ax.annotate(txt, (di.T[0][i],di.T[1][i]))
		plt.scatter(di.T[0], di.T[1], c='k', s=20)
		plt.savefig(self.outdir+'RAMA_'+name.split('/')[1] + self.name_mod + ".png")
		return None

	#####################################
	#####      Polymer Models      ######
	#####################################

	# this really needs to be cleaned up

	def gaussian_chain_2(self):
		"""
		Compute distributions and fit them to Gaussian. If fit is good, it's a Gaussian chain

		"""
		from scipy.optimize import curve_fit
		def pdf(rL,R,a):
			return (3.0/(2*np.pi*R**2))**(3.0/2.0) * np.exp(-1.5*((rL-a)/R)**2)
		#---EED
		plt.clf()
		n, bins, patches = plt.hist(self.EED, 100)
		#nbins = []
		#for b in range(len(bins)-1):
		#	nbins.append((bins[b]+bins[b+1])/2)
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		popt, pcov = curve_fit(pdf, nbins, n)
		plt.plot(nbins, pdf(nbins, *popt), 'r-')
		plt.savefig(self.outdir+'gaussian_chain_EED' + self.name_mod + '.png')
		np.savetxt(self.outdir+'gaussian_chain_EED_fit' + self.name_mod + '.npy',[nbins,pdf(nbins, *popt)])

		#--- Rg
		plt.clf()
		n, bins, patches = plt.hist(self.Rg, 100)
		#nbins = []
		#for b in range(len(bins)-1):
		#	nbins.append((bins[b]+bins[b+1])/2)
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		popt, pcov = curve_fit(pdf, nbins, n)
		plt.plot(nbins, pdf(nbins, *popt), 'r-')
		plt.savefig(self.outdir+'gaussian_chain_Rg' + self.name_mod + '.png')
		np.savetxt(self.outdir+'gaussian_chain_Rg_fit' + self.name_mod + '.npy',[nbins,pdf(nbins, *popt)])
		return None

	def gaussian_chain(self):
		"""
		Legacy 

		"""
		import matplotlib.mlab as mlab
		import scipy.stats as sci
		nbins = 30
		n, bins, patches = plt.hist(self.EED, 30, normed=1)
		(mu, sigma) = sci.norm.fit(self.EED)
		fit = mlab.normpdf(bins, mu, sigma)
		l = plt.plot(bins, fit, 'r--', linewidth=2)
		plt.show()
		return None

	def semiflexible_chain(self):
		"""
		Compute distributions and compare to skewed Gaussian. If fit is good then it's semi-flexible (like DNA)

		"""
		# Gaussian Chain
		plt.clf()
		import matplotlib.mlab as mlab
		import scipy.stats as sci
		nbins = 30
		n, bins, patches = plt.hist(self.EED, 30, normed=1)
		(mu, sigma) = sci.norm.fit(self.EED)
		fit = mlab.normpdf(bins, mu, sigma)
		l = plt.plot(bins, fit, 'r--', linewidth=2)
		#plt.show()

		# Semiflexible Chain
		from scipy.optimize import curve_fit
		def semiflex(y,Nc,p):
			return Nc * (y)**(-1.5) * (y+1)**(-3) * np.exp(-1.5 * (p/y))
		nbins = [(bins[b]+bins[b+1])/2 for b in range(len(bins)-1)]
		nbins = np.asarray(nbins).ravel()
		popt, pcov = curve_fit(semiflex, nbins, n)
		plt.plot(nbins, semiflex(nbins, *popt), 'k--')
		plt.xlabel("EED")
		plt.ylabel("P(EED)")

		#plt.show()
		plt.savefig(self.outdir+'semiflexible_chain_EED' + self.name_mod + '.png')
		np.savetxt(self.outdir+ 'semiflexible_chain_EED_fit' + self.name_mod + '.npy',[nbins,semiflex(nbins, *popt)])
		return None

	def compute_flory(self,struc):
		"""
		The coordinates need to be centered EACH TIME Rg is computed. Therefore you can't just keep adding to S and 
		recomputing the eigenvalues without centering everything. 

		"""

		N = range(5,25)
		rg = np.zeros(len(N))
		count = np.zeros(len(N))
		for n in N:
			for r in range(self.nres-n):
				sel = struc.atom_slice(struc.topology.select('resid ' + str(r) + ' to ' + str(r+n-1)))
				rg[n-5] += md.compute_rg(sel)[0] 
				count[n-5] += 1
		rg = [rg[i]/count[i] for i in range(len(rg))]
		self.fex.append(rg)
		return None

	def flory_alt(self,struc):
		"""
		The coordinates need to be centered EACH TIME Rg is computed. Therefore you can't just keep adding to S and 
		recomputing the eigenvalues without centering everything. 

		"""

		#coors = struc.xyz[0] 
		CA = struc.atom_slice(struc.topology.select('name CA'))[0]
		coors = CA.xyz[0]
		#import ipdb; ipdb.set_trace()

		def RG(coors):
        	        coors -= np.mean(coors,axis=0)
        	        S = np.sum([np.einsum('i,j->ij',c,c) for c in coors], axis=0)/(len(coors))
        	        L = np.linalg.eig(S)[0]
			return np.sqrt(np.sum(L))

		nres = len(coors)
		smallest_pep = 5
		largest_pep = 25
		N = range(smallest_pep,largest_pep)
		rg = np.zeros(len(N))
		count = np.zeros(len(N))

		for res in range(nres):
			save_coors = [rescoor for rescoor in coors[res:res+N[0]]]
			for peplen in N:
				end = res+peplen
				if end < nres:
					rescoor = coors[end]
					save_coors.append(rescoor)
					rg[peplen-smallest_pep] += RG(save_coors)
					count[peplen-smallest_pep] += 1

		print "rg2", rg
		rg = [rg[i]/count[i] for i in range(len(rg))]
		self.fex.append(rg)
		return None

	def fex_hist(self):
		"""
		Plot simple histogram of Flory exponents

		"""
		from scipy.optimize import curve_fit
		np.savetxt(self.outdir+ 'flory_exponents' + self.name_mod + '.npy',self.fex)
		def f(N,b,v):
			return N*v + b

		nframes, npep = np.array(self.fex).shape
		N = np.log(np.asarray(range(5,25)).ravel())
		rg = np.log(np.sum(self.fex,axis=0)/nframes)

		plt.clf()
		plt.plot(N,rg,color='b')
		popt, pcov = curve_fit(f, N, rg)
		print popt

		err = np.zeros(len(N)) 
		for i in range(len(N)):
			err[i] = np.std(rg[i])
		plt.errorbar(N, rg,yerr=err/(np.sqrt(nframes)-1), fmt='o',color='b')

		plt.plot(N, f(N, *popt), 'r-')
		plt.xlabel('log(N)')
		plt.ylabel('log(Rg)')
		plt.text(1.6,-0.14, r'$R_g$ = %.2f*$N^{%.2f}$' % (np.exp(popt[0]),popt[1]), fontdict=font)
		plt.subplots_adjust(top=0.9,bottom=0.15,left=0.18,right=0.85,hspace=0.2,wspace=0.2)
		plt.savefig(self.outdir+ 'flory_exponents' + self.name_mod + '.png')
		return None

	#####################################
	#####      Task Managers       ######
	#####################################
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

	def check_input(self):
		"""
		This function determines which calcs need to be run

		"""
		# make sure 'outdir' is formatted properly
		try:
			# sometimes I accidentally put / on the end twice
			if self.outdir[-2] == '/' and self.outdir[-1] == '/':
				self.outdir = self.outdir[:-1]
			# sometimes I over-correct and don't put any
			if self.outdir[-1] != '/':
				self.outdir = self.outdir+'/'
		except:
			if self.outdir == '':
				pass
			else:
				print "something is weird here. outdir =", self.outdir

		# if no calculations are specified, run these
		if self.calcs == []:
			for c in ['Gyr', 'Rg', 'SASA', 'EED', 'Asph', 'cmaps', 'SS', 'PCA']:
				if c not in self.calcs:
					self.calcs.append(c)

		# Rg and Asph come from gyration tensor
		if 'Rg' in self.calcs or 'Asph' in self.calcs:
			for c in ['Gyr']:
				if c not in self.calcs:
					self.calcs.append(c)

		if 'surface_contacts' in self.calcs:
			for c in ['SASA']:
				if c not in self.calcs:
					self.calcs.append(c)

		# PCA requires the following calculations. Make sure they're there
		elif 'PCA' in self.calcs:
			for c in ['Gyr', 'Rg', 'SASA', 'EED', 'Asph']:
				if c not in self.calcs:
					self.calcs.append(c)

		elif 'chain' in self.calcs:
			for c in ['Gyr', 'Rg', 'EED']:
				if c not in self.calcs:
					self.calcs.append(c)

		# make sure all calculations supplied are valid
		for c in self.calcs:
			if c not in ['Rg', 'SASA', 'EED', 'Asph', 'rama', 'cmaps', 'PCA', 'gcmaps','XRg','SS', 'chain', 'score','flory', 'centroids', 'Gyr', 'surface_contacts']:
				print c, "is not a known calculation. Exiting..."
				exit()

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
		Load pre-computed data from previous runs

		"""
		self.calcs = np.array(self.calcs)

		if 'Rg' in self.calcs and os.path.isfile(self.outdir+'Rg'+self.name_mod+file_ext):
			print 'loading data for Rg...', self.outdir+'Rg'+self.name_mod+file_ext
			self.Rg = np.loadtxt(self.outdir+'Rg'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'Rg')]
		if 'EED' in self.calcs and os.path.isfile(self.outdir+'EED'+self.name_mod+file_ext):
			print 'loading data for EED...', self.outdir+'EED'+self.name_mod+file_ext
			self.EED = np.loadtxt(self.outdir+'EED'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'EED')]
		if 'Asph' in self.calcs and os.path.isfile(self.outdir+'Asph'+self.name_mod+file_ext):
			print 'loading data for Asph...', self.outdir+'Asph'+self.name_mod+file_ext
			self.Asph = np.loadtxt(self.outdir+'Asph'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'Asph')]
		if 'SASA' in self.calcs and os.path.isfile(self.outdir+'SASA'+self.name_mod+file_ext):
			print 'loading data for SASA...', self.outdir+'SASA'+self.name_mod+file_ext
			self.SASA = np.loadtxt(self.outdir+'SASA'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'SASA')]
		if 'cmaps' in self.calcs and os.path.isfile(self.outdir+'CMAPS'+self.name_mod+file_ext):
			print 'loading data for CMAPS...', self.outdir+'CMAPS'+self.name_mod+file_ext
			cmaps_raw = np.loadtxt(self.outdir+'CMAPS'+self.name_mod+file_ext)
			nres = np.sqrt(cmaps_raw.shape[1]).astype(int)
			self.cmaps = cmaps_raw.reshape(cmaps_raw.shape[0], nres, nres)
			self.calcs = self.calcs[np.where(self.calcs != 'cmaps')]
		if 'gcmaps' in self.calcs and os.path.isfile(self.outdir+'GCMAPS'+self.name_mod+file_ext):
			print 'loading data for GCMAPS...', self.outdir+'GCMAPS'+self.name_mod+file_ext
			self.gcmaps = np.loadtxt(self.outdir+'GCMAPS'+self.name_mod+file_ext)
			self.calcs = self.calcs[np.where(self.calcs != 'gcmaps')]
		if 'rama' in self.calcs and os.path.isfile(self.outdir+"RAMA_all" + self.name_mod + ".npy"):
			print 'loading dihedrals...', self.outdir+"RAMA_all" + self.name_mod + ".npy"
			self.dihedrals = np.loadtxt(self.outdir+"RAMA_all" + self.name_mod + ".npy")
			self.calcs = self.calcs[np.where(self.calcs != 'rama')]
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
				print "loading PCA data..."
				print pca1
				print pca2
				print pca3
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
		if 'Rg' in self.calcs:
			np.savetxt(self.outdir+'Rg'+self.name_mod+file_ext,self.Rg)
		if 'EED' in self.calcs:
			np.savetxt(self.outdir+'EED'+self.name_mod+file_ext,self.EED)
		if 'Asph' in self.calcs:
			np.savetxt(self.outdir+'Asph'+self.name_mod+file_ext,self.Asph)
		if 'SASA' in self.calcs:
			np.savetxt(self.outdir+'SASA'+self.name_mod+file_ext,self.SASA)
		# cmaps and SS require more processing before they can be saved. They will be saved in the av_cmaps & av_SS functions
		return None

	def calculations(self,struc):
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
			#self.Rg.append(md.compute_rg(struc)[0])
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
		if 'gcmaps' in self.calcs:
			self.gremlin_contact_maps(dist)
		if 'SS' in self.calcs:
			self.SS.append(md.compute_dssp(struc))
		if 'flory' in self.calcs:
			self.compute_flory(struc)
		if 'rama' in self.calcs:
			self.compute_phipsi(struc)
		if 'surface_contacts' in self.calcs:
			self.resnames = [struc.atom_slice(struc.topology.select('name CA')).topology.atom(r).residue.name for r in range(self.nres)]
			self.surface_contacts(struc,SASA)
		return None

	def run(self,mode='default'):
	        """
		Runs and handles all function calls. All data is stored in class object.

	        """
		from timeit import default_timer as timer

		start = timer()

		global file_ext ; file_ext = '_raw.npy'

		self.mode = mode

		#---Check to see which calculations need to be run
		self.check_input()

		print self.trajname

		#---Load existing data
		if self.mode == 'default' and self.OVERWRITE == 'n': 
			self.load_data()
		elif self.OVERWRITE == 'y':
			print "OVERWRITING OLD DATA!"

		print "calculations left to do:", self.calcs

		#---Decide if it's necessary to load trajectories/PDBs
		if len(self.calcs) == 0: LOAD = False
		else: LOAD = True

		if self.mode == 'default':
			if LOAD == True:
				print "Loading Trajectory"
				#---Run Calculations *DCD*
				if self.trajname.split('.')[-1] == 'dcd':
					dcd = md.load_dcd(self.trajname, self.top)
					for struc in dcd:
						self.calculations(struc)

				#                    *PDB*
				elif self.trajname.split('.')[-1] == 'txt':
					for PDB in open(self.trajname):
						struc = md.load(PDB.split("\n")[0])
						self.calculations(struc)
					self.top = PDB.split("\n")[0]

				#---Write data
				self.write_data()

		#---Code can be run in special mode where it references a Rosetta score file and only computes statistics on top N structures
		elif self.mode == 'score':
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
				scores = self.load_scores(self.score_file) 
			#---Save all Rosetta ID's with their associated scores
			frame_sel = []
			ind = np.argsort(scores)[:self.ros_frames]
			for i in ind:
				self.scores.append([i, scores[i]])
				frame_sel.append(i)
			#---Save associated files in a separate directory
			fr = 0
			if not os.path.exists(self.outdir+"top_score/pdbs"):
				os.makedirs(self.outdir+"top_score/pdbs")
			#---Iterate through pdb's with best scores
			for PDB in open(self.trajname):
				if fr in frame_sel:
					struc = md.load(PDB.split("\n")[0])
					struc.save(self.outdir+"top_score/"+PDB.split("\n")[0].split("/")[-1])
					self.calculations(struc)
					print PDB.split("\n")[0]
				fr += 1
			self.top = PDB.split("\n")[0]

		#---Post-Processing
		if 'cmaps' in self.calcs:
			try: self.av_cmaps(self.cmaps)
			except: print "CMAPS didnt work"
		if 'gcmaps' in self.calcs:
			try: self.av_cmaps(self.gcmaps,"gremlin")
			except: print "grem CMAPS didnt work"
		if 'surface_contacts' in self.calcs:
			self.av_cmaps(self.scmaps,"surface")
			#try: self.av_cmaps(self.scmaps,"surface")
			#except: print "surface CMAPS didnt work"
		if 'SS' in self.calcs:
			try: self.av_SS()
			except: print "SS didnt work" ; exit()
		if 'EED' in self.calcs and 'Asph' in self.calcs:
			try: self.scatterplot(self.EED, self.Asph, 'EED', 'Asph', 'EED_v_Asph')
			except: print "didnt work 3"
		if 'Rg' in self.calcs and 'SASA' in self.calcs:
			try: self.scatterplot(self.Rg, self.SASA, 'Rg', 'SASA', 'Rg_v_SASA')
			except: print "didnt work 4"
		if 'PCA' in self.calcs:
			self.run_PCA()
		if 'flory' in self.calcs:
			self.fex_hist()
		if 'chain' in self.calcs:
			#self.gaussian_chain()
			self.semiflexible_chain()
		if 'centroids' in self.calcs:
			self.cluster_centroids()
		if 'rama' in self.calcs:
			self.rama()

		end = timer()
		print "Total execution time:", end-start

		return None

def USAGE():
	print "USEAGE: python rosetta_analysis.py ARGS"
	print "ARGS: EITHER a list of pdbs in a file with a .txt extension, or"
	print "a .dcd and a .pdb"
	exit()

if __name__ == "__main__":
	import sys
	if len(sys.argv) == 2:
		if sys.argv[1].split('.')[1] == 'txt':
			PDB_list = sys.argv[1]
			sa = SA(PDB_list,'','_test','',['surface_contacts'])
		else:
			USAGE()
	elif len(sys.argv) == 3:
		if sys.argv[1].split('.')[1] == 'dcd' and sys.argv[2].split('.')[1] == 'pdb':
			dcd = sys.argv[1]
			pdb = sys.argv[2]
			sa = SA(dcd,pdb,'test','test_traj/',[])
		else:
			USAGE()
	else:
		USAGE()

	sa.overwrite()
	sa.run('score')
	#sa.ros_score_sort()
	#sa.write_clusters()
	#sa.analyze_clusters()
