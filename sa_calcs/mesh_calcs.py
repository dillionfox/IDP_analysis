from utils import np, md, plt

scale = 2
class mesh_calcs:
	"""
	Lipids are forced onto a 2D rectangular mesh in the XY plane

	Calcs:
	- Heights (av_heights)
		Z-heights of lipids are averaged over time. I added an extra step to
		apply light Gaussian smoothing but this is only for aesthetics. 
	- Thickness (thickness)
		Average thickness of bilayer. It is assumed that leaflets are
		symmetric and lipid IDs with lower numbers or on upper leaflet
	- Interdigitation
		Distances between distal carbons in acyl chains are averaged.
		This is one way to detect interdigitated gel phase.

	I hard-coded some things in here, so be careful. It shouldn't be hard
	to change them.

	"""

	def __init__(self,outdir,name_mod):
		self.area_per_lipid = [] # Area Per Lipid
		self.outdir = outdir
		self.av_heights = np.array([])
		self.n_box_counts = np.array([])
		self.av_prot_coors = []
		self.lipid_coors = {}
		self.lipid_ids = []
		self.n_lipids = None
		self.nx = None # dimensions of box
		self.ny = None # dimensions of box 
		self.cx = None # center of box
		self.cy = None # center of box 
		self.av_di = None

	@staticmethod
	def snap_to_grid(q):
		"""
		Convert xyz coordinates to matrix indices

		"""
		return int(np.round(scale*q))-1

	def wrap_back(self,c,d):
		"""
		Keep all coordinates in pbcs

		"""
		if c >= d:
			return c-d
		elif c < 0:
			return c+d
		return c

	def get_shift(self,fr):
		"""
		Shift each frame so the protein is in the center

		"""
		prot_coor = self.av_prot_coors[fr]
		center_x, center_y = self.snap_to_grid(prot_coor[0]),self.snap_to_grid(prot_coor[1])
		return center_x-self.cx,center_y-self.cy

	# core functions
	def interdigit(self,struc,fr):
		"""
		Average interdigitation of acyl tails

		"""
		shift_x, shift_y = self.get_shift(fr)
		for ind,lipid in enumerate(self.lipid_ids):
			lc = self.lipid_coors['interdigitation'][ind][fr]
			x = self.wrap_back(self.snap_to_grid(lc[0])-shift_x,self.nx) 
			y = self.wrap_back(self.snap_to_grid(lc[1])-shift_y,self.ny)
			z = float(lc[2])
			if ind < self.n_lipids//2:
				leaflet = 0 # upper leaflet
			else:
				leaflet = 1 # lower leaflet
			self.interdigitation[leaflet][fr][x][y] += z
			self.interdigit_count[leaflet][fr][x][y] += 1
		return None

	def heights(self,struc,fr):
		"""
		Average height of BOTH leaflets

		"""
		shift_x, shift_y = self.get_shift(fr)
		for ind,lipid in enumerate(self.lipid_ids):
			lc = self.lipid_coors['heights'][ind][fr]
			x = self.wrap_back(self.snap_to_grid(lc[0])-shift_x,self.nx) 
			y = self.wrap_back(self.snap_to_grid(lc[1])-shift_y,self.ny)
			z = np.mean(lc[2])
			self.n_box_counts[x][y] += 1
			self.av_heights[x][y] += z
		return None

	def thickness(self,struc,fr):
		"""
		Average thickness of the membrane as f(x,y)

		"""
		shift_x, shift_y = self.get_shift(fr)
		for ind,lipid in enumerate(self.lipid_ids):
			lc = self.lipid_coors['thickness'][ind][fr]
			x = self.wrap_back(self.snap_to_grid(lc[0])-shift_x,self.nx) 
			y = self.wrap_back(self.snap_to_grid(lc[1])-shift_y,self.ny)
			z = np.mean(lc[2])
			if ind < self.n_lipids//2:
				leaflet = 0 # upper leaflet
			else:
				leaflet = 1 # lower leaflet
			self.d_z[leaflet][fr][x][y] += z
			self.dz_count[leaflet][fr][x][y] += 1
		return None

	# post-analysis functions
	def plot_interdigitation(self):
		"""
		Plot thickness of membrane and try to apply smoothing (though it is usually not necessary)

		"""
		av_di_u = np.divide(np.sum(self.interdigitation[0],axis=0),np.sum(self.interdigit_count[0],axis=0))
		av_di_l = np.divide(np.sum(self.interdigitation[1],axis=0),np.sum(self.interdigit_count[1],axis=0))
		self.av_di = av_di_u - av_di_l
		np.savetxt(self.outdir+'average_interdigitation_centered.npy',self.av_di)
		fig, ax = plt.subplots()
		cax = ax.imshow(self.av_di,vmin=-1.5,vmax=0.5)
		fig.colorbar(cax)
		ax = self.add_prot_frame(ax)
		print "plotting", self.outdir+'average_interdigitation_centered.png'
		plt.savefig(self.outdir+'average_interdigitation_centered.png')
		try:
			self.artificial_smoothing(self.interdigitation,'average_interdigitation_centered_smoothed_s-')
		except:
			print "artificial smoothing didn't work"
		return None

	def normalize_heights(self):
		"""
		Average heights over time

		"""
		for r,row in enumerate(self.av_heights):
			for c in range(len(row)):
				self.av_heights[r][c]/=float(self.n_box_counts[r][c])
		fig, ax = plt.subplots()
		cax = ax.imshow(self.av_heights)
		fig.colorbar(cax)
		print "plotting", self.outdir+'average_heights_centered.png'
		plt.savefig(self.outdir+'average_heights_centered.png')
		self.artificial_smoothing(self.av_heights,'average_heights_centered_smoothed_s-')
		return None

	def plot_thickness(self):
		"""
		Plot thickness of membrane and try to apply smoothing (though it is usually not necessary)

		"""
		self.av_dz_u = np.divide(np.sum(self.d_z[0],axis=0),np.sum(self.dz_count[0],axis=0))
		self.av_dz_l = np.divide(np.sum(self.d_z[1],axis=0),np.sum(self.dz_count[1],axis=0))
		self.av_dz = self.av_dz_u - self.av_dz_l
		fig, ax = plt.subplots()
		cax = ax.imshow(self.av_dz)
		fig.colorbar(cax)
		print "plotting", self.outdir+'average_thickness_centered.png'
		plt.savefig(self.outdir+'average_thickness_centered.png')
		try:
			self.artificial_smoothing(self.thickness,'average_heights_centered_smoothed_s-')
		except:
			print "artificial smoothing didn't work"
		return None

	def add_prot_frame(self,ax):
		"""
		Add outline of protein from last frame to plot

		"""
		prot_coor = self.CA_coors[-1]
		COM = np.mean(prot_coor,axis=0)
		com_coors = prot_coor-COM+[self.nx/2.0,self.ny/2.0,0]
		ax.scatter(com_coors.T[0],com_coors.T[1],facecolor='None',edgecolor='w')
		return ax

	@staticmethod
	def artificial_smoothing(im,name):
		import scipy
		import scipy.ndimage
		for s in [1]:
			plt.clf()
			f0 = scipy.ndimage.gaussian_filter(im,sigma=s,order=0,mode='wrap')
			fig, ax = plt.subplots()
			cax = ax.imshow(f0)
			fig.colorbar(cax)
			ax.scatter(px, py, s=500, marker='s',facecolor='None',edgecolor='w')
			plt.savefig(self.outdir+name+str(s)+'.png')
		return None

	# pre-calcs
	def mesh_precalcs(self,traj):
		self.n_frames = len(traj)
		struc = traj[0]
		sel_ind = struc.topology.select('resname DMPC and name C21')	# select one atom per lipid to extract residue IDs
		sel = struc.atom_slice(sel_ind)					# MDTraj object to store trajectory information
		table, bonds = sel.topology.to_dataframe()			# Extract info from MDTraj object
		self.lipid_ids = [i for i in table['resSeq']]			# residue IDs of each lipid
		self.n_lipids = len(self.lipid_ids)				# total number of lipids
		self.nx = scale*int(np.ceil(struc.unitcell_lengths[0][0]))	# number of indices in x
		self.ny = scale*int(np.ceil(struc.unitcell_lengths[0][1]))	# number of indices in y
		self.av_heights = np.zeros((self.nx,self.ny))			# average height of bilayer midplane
		self.n_box_counts = np.zeros((self.nx,self.ny))			# empty matrix to keep track of av_heights counts
		self.d_z = np.zeros((2,self.n_frames,self.nx,self.ny))		# thickness of membrane
		self.dz_count = np.zeros((2,self.n_frames,self.nx,self.ny))	# thickness of membrane
		self.interdigitation = np.zeros((2,self.n_frames,self.nx,self.ny))	# thickness of membrane
		self.interdigit_count = np.zeros((2,self.n_frames,self.nx,self.ny))	# thickness of membrane
		self.av_dz_u = np.zeros((self.nx,self.ny))			# thickness of membrane
		self.av_dz_l = np.zeros((self.nx,self.ny))			# thickness of membrane
		self.cx = self.nx/2						# center of box, x (in discrete mesh coordinates)
		self.cy = self.ny/2						# center of box, y (in discrete mesh coordinates) 
		return None

	def protein_mesh(self,traj):
		self.protein_coors = traj.atom_slice(traj.topology.select('protein and backbone')).xyz
		self.av_prot_coors = np.mean(self.protein_coors,axis=1)
		self.CA_coors = traj.atom_slice(traj.topology.select('protein and name CA')).xyz
		return None

	def lipid_mesh(self,traj,calc,hg_atoms):
		"""
		Initiate mesh generation by extracting coordinates

		"""
		# lipid COMs
		coors = []
		for lipid in self.lipid_ids:
			lc = np.mean(traj.atom_slice(traj.topology.select(hg_atoms+' and residue ' +str(lipid))).xyz,axis=1)
			coors.append(lc)
		self.lipid_coors[calc] = np.array(coors)

		# protein COM
		return None
