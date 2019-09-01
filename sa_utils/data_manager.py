import numpy as np
import os

file_ext = '_raw.npy'

class data_manager:

	def __init__(self,sa_dict):
		self.__dict__.update((key, sa_dict[key]) for key in sa_dict)

	def load_data(self):
		"""
		Don't compute things twice. Load pre-computed data from previous runs

		"""
		#---Load data
		for c in ['Rg', 'EED', 'Asph', 'SASA', 'cmaps', 'gcmaps', 'rama', 'MASA', 'diffusion', 'flory', 'rmsd',\
				'membrane_contacts','av_heights','av_interdigitation']:
			if c in self.calcs_master and os.path.isfile(self.outdir+c+self.name_mod+file_ext):
				print 'loading data for', c, self.outdir+c+self.name_mod+file_ext
				self.__dict__[c] = np.loadtxt(self.outdir+c+self.name_mod+file_ext)
				self.calcs = self.calcs[np.where(self.calcs != c)]
				self.traj_list = self.traj_list[np.where(self.traj_list != c)]

		if 'membrane_contacts' in self.calcs_master and os.path.isfile(self.outdir+'contact_frames'+self.name_mod+file_ext):
			self.contact_frames = np.loadtxt(self.outdir+'contact_frames'+self.name_mod+file_ext)

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
		for c in ['Rg', 'EED', 'Asph', 'SASA', 'rama', 'diffusion', 'flory', 'rmsd', 'membrane_contacts',\
				'area_per_lipid','av_heights','av_interdigitation']:
			if c in self.calcs_master:
				print "Writing data for", c, ":", self.outdir+c+self.name_mod+file_ext
				try:
					np.savetxt(self.outdir+c+self.name_mod+file_ext,self.__dict__[c])
				except:
					import ipdb; ipdb.set_trace()
		if 'membrane_contacts' in self.calcs_master:
			np.savetxt(self.outdir+'contact_frames'+self.name_mod+file_ext,self.contact_frames)
		#if 'interdigitation' in self.calcs_master:
		#	np.savetxt(self.outdir+c+self.name_mod+file_ext,self.__dict__['av_di'])
		return None

