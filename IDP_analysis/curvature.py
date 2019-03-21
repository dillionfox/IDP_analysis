from lib_handler import np, md, plt

scale = 2
hg_atoms = 'resname DMPC and not type H'
class curvature:
	"""
	Right now the name "curvature" is more of a goal than a description.
	This class computes the average heights of lipids in a bilayer. 

	The lipids are forced onto a 2D rectangular mesh in the XY plane
	and the heights are averaged over time. I added an extra step to
	apply light Gaussian smoothing but this is only for aesthetics. 

	I hard-coded some things in here, so be careful. It shouldn't be hard
	to change them.

	"""

	def __init__(self,outdir,name_mod):
		self.area_per_lipid = [] # Area Per Lipid
		self.outdir = outdir
		self.av_heights = np.array([])
		self.n_box_counts = np.array([])
		self.protein_coors = []
		self.lipid_coors = []
		self.lipid_ids = []
		self.n_lipids = None
		self.nx = None # dimensions of box
		self.ny = None # dimensions of box 
		self.cx = None # center of box
		self.cy = None # center of box 

	@staticmethod
	def snap_to_grid(q):
		return int(np.round(scale*q))-1

	def wrap_back(self,c,d):
		if c >= d:
			return c-d
		elif c < 0:
			return c+d
		return c

	def heights(self,struc,fr):
		prot_coor = self.protein_coors[fr]
		center_x, center_y = self.snap_to_grid(prot_coor[0]),self.snap_to_grid(prot_coor[1])
		shift_x, shift_y = center_x-self.cx,center_y-self.cy
		for ind,lipid in enumerate(self.lipid_ids):
			lc = self.lipid_coors[ind][fr]
			x = self.wrap_back(self.snap_to_grid(lc[0])-shift_x,self.nx) 
			y = self.wrap_back(self.snap_to_grid(lc[1])-shift_y,self.ny)
			z = np.mean(lc[2])
			self.n_box_counts[x][y] += 1
			self.av_heights[x][y] += z
		return None

	def normalize_heights(self):
		for r,row in enumerate(self.av_heights):
			for c in range(len(row)):
				self.av_heights[r][c]/=float(self.n_box_counts[r][c])
		prot_coor = self.protein_coors[0][0:2]
		center_x, center_y = self.snap_to_grid(prot_coor[0]),self.snap_to_grid(prot_coor[1])
		shift_x, shift_y = center_x-self.cx,center_y-self.cy
		px,py = center_x-shift_x,center_y-shift_y
		fig, ax = plt.subplots()
		cax = ax.imshow(self.av_heights)
		fig.colorbar(cax)
		ax.scatter(px, py, s=500, marker='s',facecolor='None',edgecolor='w')
		print "plotting", self.outdir+'average_heights_centered.png'
		plt.savefig(self.outdir+'average_heights_centered.png')

		import scipy
		import scipy.ndimage
		for s in [1]:
			plt.clf()
			av_heights_f0 = scipy.ndimage.gaussian_filter(self.av_heights,sigma=s,order=0,mode='wrap')
			fig, ax = plt.subplots()
			cax = ax.imshow(av_heights_f0)
			fig.colorbar(cax)
			ax.scatter(px, py, s=500, marker='s',facecolor='None',edgecolor='w')
			plt.savefig(self.outdir+'average_heights_centered_smoothed_s-'+str(s)+'.png')
		return None

	def heights_precalcs(self,struc):
		sel_ind = struc.topology.select('resname DMPC and name C21')
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		self.lipid_ids = [i for i in table['resSeq']]
		self.n_lipids = len(self.lipid_ids)
		#self.lipid_ids = self.lipid_ids[:self.n_lipids/2] # upper leaflet only
		self.nx = scale*int(np.ceil(struc.unitcell_lengths[0][0]))
		self.ny = scale*int(np.ceil(struc.unitcell_lengths[0][1]))
		self.av_heights = np.zeros((self.nx,self.ny))
		self.n_box_counts = np.zeros((self.nx,self.ny))
		self.cx = self.nx/2
		self.cy = self.ny/2
		return None

	def lipid_mesh(self,traj):
		# lipid COMs
		for lipid in self.lipid_ids:
			lc = np.mean(traj.atom_slice(traj.topology.select(hg_atoms+' and residue ' +str(lipid))).xyz,axis=1)
			self.lipid_coors.append(lc)
		self.lipid_coors = np.array(self.lipid_coors)

		# protein COM
		self.protein_coors = np.mean(traj.atom_slice(traj.topology.select('protein and backbone')).xyz,axis=1)
		return None
