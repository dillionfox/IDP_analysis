from utils import np, md, os

#from sa import md, os, np

def write_clusters(k,trajname,top,outdir,name_mod):
	"""
	This function writes a new dcd trajectory for each k-means cluster defined by PCA functions.
	This might be useful if you want to perform structural analyses on individual clusters.

	"""
	print 'loaded', outdir+'PCA_kmeans_clusters' + name_mod + '.npy'
	k = np.loadtxt(outdir+'PCA_kmeans_clusters' + name_mod + '.npy')

	# Check to see if cluster trajectories already exist
	count = 0
	for cluster in range(int(max(k))+1):
		if os.path.isfile(outdir + 'cluster_' + str(cluster) + name_mod + '.dcd'):
			count+=1

	if count == int(max(k)+1):
		print "cluster trajectories already exist"
		return None

	#---DCD
	if trajname.split('.')[-1] == 'dcd':
		dcd = md.load_dcd(trajname, top)
		for cluster in range(int(max(k))+1):
			with md.formats.DCDTrajectoryFile(outdir + 'cluster_' + str(cluster) + name_mod + '.dcd', 'w') as f: 
				for frames in np.where(k==cluster):
					for fr in frames:
						f.write(dcd[fr].xyz*10)

	#---PDB
	elif trajname.split('.')[-1] == 'txt':
		struc = []
		for PDB in open(trajname):
			struc.append(md.load(PDB.split("\n")[0]))

		for cluster in range(int(max(k))+1):
			with md.formats.DCDTrajectoryFile(outdir + 'cluster_' + str(cluster) + name_mod + '.dcd', 'w') as f: 
				for frames in np.where(k==cluster):
					for fr in frames:
						f.write(struc[fr].xyz*10)
				print "done with cluster", cluster
	return None

