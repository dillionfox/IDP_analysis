from lib_handler import np, plt

"""
run_PCA(EED,Rg,SASA,Asph,outdir,name_mod,mode,scores,trajname,ros_frames): None
compute_PCA(A): P
plot_Kmeans_PCA(A, B, ind_list, PC_labels, kmeansoutput, outdir, name_mod, scores, ros_frames): None
extract_PCA(pca, PC_labels, mode, outdir, name_mod): ind_list
compute_Kmeans(B, mode): j
norm_shift(vec): vec
extract_inds(fil): PC_labels 1 & 2
cluster_centroids(B,k,calcs,mode,scores,outdir,name_mod,trajname): None

"""

class pca:
	def __init__(self):
		self.k = []
		self.B = []
	
	def plot_Kmeans_PCA(self,ind_list, PC_labels, kmeansoutput, outdir, name_mod, scores, ros_frames):
		"""
		Make plot in PC space with points colored by k-means
	
		"""
		plt.clf()
		plt.figure('Top 2 PCAs colored by K-Means Clustering')
		compare_2 = 'n'
		#---Plot in PCA Space
		if compare_2 == 'y':
			plt.scatter(self.B[:,0][:50], self.B[:,1][:50], c=kmeansoutput.labels_[:50],edgecolor='k')
			plt.scatter(self.B[:,0][50:], self.B[:,1][50:], c=kmeansoutput.labels_[50:],marker="s")
		else:
			plt.scatter(self.B[:,0], self.B[:,1], c=kmeansoutput.labels_)
	
		plt.xlabel('PC-1')
		plt.ylabel('PC-2')
		plt.xlim([-1,1])
		plt.ylim([-1,0.98])
		plt.title('Top 2 PCAs colored by K-Means Clustering')
		plt.subplots_adjust(top=0.9,bottom=0.15,left=0.18,right=0.85,hspace=0.2,wspace=0.2)
		plt.savefig(outdir+'PCA_kmeans' + name_mod + '.png')
		if name_mod == '_ros_score':
			ind = np.array(scores).T[0]
			global_labels = [int(ind[l]) for l in range(ros_frames)]
			np.savetxt(outdir+'cluster_labels'+name_mod+ '.txt',global_labels,fmt="%s")
		np.savetxt(outdir+'PCA_kmeans_clusters' + name_mod + '.npy',np.array(kmeansoutput.labels_))
		np.savetxt(outdir+'PCA_kmeans_B' + name_mod + '.npy',np.array(self.B))
		return None
	
	@staticmethod
	def extract_PCA(pca, PC_labels, mode, outdir, name_mod):
		"""
		Figure out which components are the most important
	
		"""
		elbow_cutoff = 0.85
		if mode == 'score':
			elbow_cutoff = 0.8
	
		print "PCA Variance:", pca.explained_variance_ratio_.cumsum()
		i = 1
		for pc in pca.explained_variance_ratio_.cumsum():
			if float(pc) > elbow_cutoff and i >1:
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
	
		fil = open(outdir+'PCA'+name_mod+ '.txt','w')
		fil.write("Important Components:")
		for label in range(len(ind_list)):
			fil.write(PC_labels[ind_list[label]]+",")
	
		return ind_list
	
	def compute_Kmeans(self,mode):
		"""
		Cluster the points in PC space
	
		"""
		from sklearn.cluster import KMeans
		elbow_cutoff = 0.7
		if mode == 'score':
			elbow_cutoff = 0.75
		Nc = range(1, 20)
		kmeans = [KMeans(n_clusters=i) for i in Nc]
		score = [kmeans[i].fit(self.B).score(self.B) for i in range(len(kmeans))]
	
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
	
	@staticmethod
	def norm_shift(vec):
		"""
		Force all data to extend from -1 to 1
	
		"""
		vec = np.array(vec)		# [a, b]
		vec -= min(vec)			# [0, b-a]
		vec /= (max(vec) - min(vec)) 	# [0, 1]
		vec *= 2			# [0, 2]
		vec -= 1			# [-1, 1]
		return vec
	
	@staticmethod
	def compute_PCA(A):
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
		#M_PCA = compute_PCA(A)
		#print "mine:\n",  np.real(M_PCA).T
		#print "sklearn\n", B
		return P
	
	def run_PCA(self,EED,Rg,SASA,Asph,outdir,name_mod,mode,scores,trajname,ros_frames,calcs):
		"""
		Main function for running PCA. This function also uses k-means clustering on the PC's.
	
		"""
		from sklearn.decomposition import PCA
		from sklearn.cluster import KMeans
	
		# normalize and shift all vectors to be centered around zero
		norm_EED =  self.norm_shift(EED)
		norm_Rg =   self.norm_shift(Rg)
		norm_SASA = self.norm_shift(SASA) 
		norm_Asph = self.norm_shift(Asph)
	
		#---Prepare array containing all SCALED data
		A = np.array([np.array(norm_EED), np.array(norm_Rg), np.array(norm_SASA), np.array(norm_Asph)]).T
		PC_labels = ['End-to-End Distance', 'Radius of Gyration', 'SASA', 'Asphericity']
	
		#---Do PCA 
		pca = PCA(len(PC_labels))
		pca.fit(A)
		self.B = pca.transform(A)
	
		#---ind_list contains the important components
		ind_list = self.extract_PCA(pca, PC_labels, mode, outdir, name_mod)
	
		#---Do initial K-means clustering to determine number of clusters
		nK = self.compute_Kmeans(mode)
	
		#---Use optimum number of clusters for k-means
		kmeans=KMeans(n_clusters=nK)
		if len(ind_list) > 1:
			kmeansoutput=kmeans.fit(np.array([self.B[:,0],self.B[:,1]]).T)
	
			#---Plot top 2 PCA clusters colored by kmeans
			self.plot_Kmeans_PCA(ind_list, PC_labels, kmeansoutput, outdir, name_mod, scores, ros_frames)
			self.k = np.array(kmeansoutput.labels_)
			self.cluster_centroids(calcs,mode,scores,outdir,name_mod,trajname) 
	
			return None
		else:
			#---Don't bother continuing if there's only one PC. Not much to do in this case...
			print "Only one PC:", PC_labels[ind_list[0]]
			print "Not generating plot, but will generate dummy files"
			np.savetxt(outdir+'PCA_kmeans_clusters' + name_mod + '.npy',np.array([]))
			np.savetxt(outdir+'PCA_kmeans_B' + name_mod + '.npy',np.array(self.B))
	
		return None
	
	@staticmethod
	def extract_inds(fil):
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
	
	def cluster_centroids(self,calcs,mode,scores,outdir,name_mod,trajname):
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
		
		if 'PCA' not in calcs:
			ind_list = self.extract_inds(outdir+'PCA'+name_mod+ '.txt')
	
		N = int(max(self.k))
		centroids = []
		labels = []
		coors = []
		self.B = np.array([self.B[:,0], self.B[:,1]]).T
		for n in range(N+1):
			b = self.B[np.where(self.k == n)]
			val = np.mean(b,axis=0)
			centroids.append(val)
			[b_i, centroid] = find_nearest(b,val)
			B_i = np.where(self.B == centroid)[0][0]
			labels.append(B_i)
			coors.append(b[b_i])
		if mode == 'score':
			ind = np.array(scores).T[0]
			for l in labels:
				print "len ind", len(ind), "l", l
			global_labels = [int(ind[l]) for l in labels]
			pdbs = [PDB.split("\n")[0] for PDB in open(trajname)]
			label_names = [pdbs[l] for l in global_labels]
			np.savetxt(outdir+'cluster_centroid_labels'+name_mod+ '.txt',label_names,fmt="%s")
		else:
			np.savetxt(outdir+'cluster_centroid_labels'+name_mod+ '.txt',labels, fmt='%i')
	
		if replot_PCA in ['y', 'yes', 'on']:
			centroids = np.array(centroids).T
			coors = np.array(coors).T
			plt.scatter(self.B[:,0], self.B[:,1], c=self.k,s=30)
			plt.scatter(centroids[0], centroids[1], c='k',s=30)
			plt.scatter(coors[0], coors[1], c='r',s=30)
			for i in range(N+1):
				plt.plot([centroids[0][i], coors[0][i]], [centroids[1][i], coors[1][i]])
			plt.xlim([-2,2])
			plt.ylim([-2,2])
			plt.xlabel('PC-1')
			plt.ylabel('PC-2')
			#plt.show()
			plt.savefig(outdir+'cluster_centroid_labels'+name_mod+ '.png')
		return None
