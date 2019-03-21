from lib_handler import np, md, plt, subprocess, os

ATOMIC_RADII = {'H'   : 0.120, 'He'  : 0.140, 'Li'  : 0.076, 'Be' : 0.059, 
                 'B'   : 0.192, 'C'   : 0.170, 'N'   : 0.155, 'O'  : 0.152, 
                 'F'   : 0.147, 'Ne'  : 0.154, 'Na'  : 0.102, 'Mg' : 0.086, 
                 'Al'  : 0.184, 'Si'  : 0.210, 'P'   : 0.180, 'S'  : 0.180, 
                 'Cl'  : 0.181, 'Ar'  : 0.188, 'K'   : 0.138, 'Ca' : 0.114, 
                 'Sc'  : 0.211, 'Ti'  : 0.200, 'V'   : 0.200, 'Cr' : 0.200, 
                 'Mn'  : 0.200, 'Fe'  : 0.200, 'Co'  : 0.200, 'Ni' : 0.163, 
                 'Cu'  : 0.140, 'Zn'  : 0.139, 'Ga'  : 0.187, 'Ge' : 0.211, 
                 'As'  : 0.185, 'Se'  : 0.190, 'Br'  : 0.185, 'Kr' : 0.202, 
                 'Rb'  : 0.303, 'Sr'  : 0.249, 'Y'   : 0.200, 'Zr' : 0.200, 
                 'Nb'  : 0.200, 'Mo'  : 0.200, 'Tc'  : 0.200, 'Ru' : 0.200, 
                 'Rh'  : 0.200, 'Pd'  : 0.163, 'Ag'  : 0.172, 'Cd' : 0.158, 
                 'In'  : 0.193, 'Sn'  : 0.217, 'Sb'  : 0.206, 'Te' : 0.206, 
                 'I'   : 0.198, 'Xe'  : 0.216, 'Cs'  : 0.167, 'Ba' : 0.149, 
                 'La'  : 0.200, 'Ce'  : 0.200, 'Pr'  : 0.200, 'Nd' : 0.200, 
                 'Pm'  : 0.200, 'Sm'  : 0.200, 'Eu'  : 0.200, 'Gd' : 0.200, 
                 'Tb'  : 0.200, 'Dy'  : 0.200, 'Ho'  : 0.200, 'Er' : 0.200, 
                 'Tm'  : 0.200, 'Yb'  : 0.200, 'Lu'  : 0.200, 'Hf' : 0.200, 
                 'Ta'  : 0.200, 'W'   : 0.200, 'Re'  : 0.200, 'Os' : 0.200, 
                 'Ir'  : 0.200, 'Pt'  : 0.175, 'Au'  : 0.166, 'Hg' : 0.155, 
                 'Tl'  : 0.196, 'Pb'  : 0.202, 'Bi'  : 0.207, 'Po' : 0.197, 
                 'At'  : 0.202, 'Rn'  : 0.220, 'Fr'  : 0.348, 'Ra' : 0.283, 
                 'Ac'  : 0.200, 'Th'  : 0.200, 'Pa'  : 0.200, 'U'  : 0.186, 
                 'Np'  : 0.200, 'Pu'  : 0.200, 'Am'  : 0.200, 'Cm' : 0.200, 
                 'Bk'  : 0.200, 'Cf'  : 0.200, 'Es'  : 0.200, 'Fm' : 0.200, 
                 'Md'  : 0.200, 'No'  : 0.200, 'Lr'  : 0.200, 'Rf' : 0.200, 
                 'Db'  : 0.200, 'Sg'  : 0.200, 'Bh'  : 0.200, 'Hs' : 0.200, 
                 'Mt'  : 0.200, 'Ds'  : 0.200, 'Rg'  : 0.200, 'Cn' : 0.200, 
                 'Uut' : 0.200, 'Fl'  : 0.200, 'Uup' : 0.200, 'Lv' : 0.200, 
                 'Uus' : 0.200, 'Uuo' : 0.200} 

class MASA:
	"""
	Membrane Accessible Surface Area

	This code can be used to very carefully identify protein-lipid contacts
	OR to calculate the surface area of the membrane.

	- MASA_protein
	Not sure if other people use this specific algorithm for other things,
	but it's based on the Shrake-Rupley SASA algorithm. Instead of searching
	for atoms that have areas exposed to no other atoms of that selection,
	this code searches for atoms of selection 1 that are exposed to atoms
	of selection 2. This is done using the same basic algorithm but the 
	conditionals are different.

	- MASA
	Shrake-Rupley

	"""

	def __init__(self,outdir):
		self.outdir = outdir
		self.sel1 = "protein"
		self.sel2 = "resname DMPC"
		#self.n_sphere_point = 960
		self.n_sphere_point = 24
		self.cutoff = 0.6
		self.probe = 0.3
		self.const = 4.0*np.pi/self.n_sphere_point
		self.sphere_points = self.generate_sphere_points()
		self.area_per_lipid = []

	def generate_sphere_points(self):
		"""
		generate coordinates on a sphere using the Golden Section Spiral algorithm
	
		"""
		points = np.zeros((self.n_sphere_point, 3))
		inc = np.pi * (3 - np.sqrt(5))
		offset = 2/float(self.n_sphere_point)
		for k in range(int(self.n_sphere_point)):
			y = k*offset - 1 + (offset/2)
			r = np.sqrt(1 - y*y)
			phi = k*inc
			points[k, :] = [np.cos(phi) * r, y, np.sin(phi) * r]
		return points
	
	@staticmethod
	def pbc_dist(r_i,r_j,pbc):
		"""
		Wrap back into box

		"""
		d = r_i - r_j
		for q in range(3):
			if abs(d[q]) > pbc[q]/2.0:
				if d[q] > 0: d[q] -= pbc[q]
				else: d[q] += pbc[q]
		return d

	@staticmethod
	def get_residues(struc,sel):
		"""
		Get residues from MDTraj object

		"""
		sel_ind = struc.topology.select(sel)
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		return table['serial'], table['resSeq'], table['resName']

	@staticmethod
	def write_pdb(coors,pdbname):
		"""
		generic pdb writing function for "test points"
	
		"""
	
		outfile = open(pdbname,"w")
		it = 0
		for i in coors:
			t1 = "ATOM"					# ATOM
			t2 = 1						# INDEX
			t3 = "C"					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = "GLY"					# RESIDUE NAME
			t6 = "A"					# CHAIN
			t7 = 0						# RESIDUE NUMBER
			t8 = ""						# INSERTION CODE
			t9 =  10*float(i[0])				# X
			t10 = 10*float(i[1])				# Y
			t11 = 10*float(i[2])				# Z
			t12 = 0.0					# OCCUPANCY
			t13 = 0.0					# TEMPERATURE FACTOR
			t14 = ""					# ELEMENT SYMBOL
			t15 = ""					# CHARGE ON ATOM
			outfile.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15))
			it+=1
		outfile.close()
		return 0

	def apl_precalcs(self,traj):
		"""
		Get some things out of the way before calculating area per lipid 
	
		"""
		struc = traj[0] ; it = 1
		sel = "resname DMPC and not type H"
		lipid_ind,lipid_res,lipid_resname = self.get_residues(struc,sel)
		lipids = traj.atom_slice(lipid_ind)

		self.r_i = []
		self.coors_i = []
		self.neighbor_indices = []

		# this atom
		for it,atom_i in enumerate(lipid_ind):
			this_atom = lipids.atom_slice(lipids.topology.select('index '+str(it)))
			self.r_i.append(self.probe + ATOMIC_RADII[this_atom.topology.to_dataframe()[0]['element'][0]])
			self.coors_i.append(this_atom.xyz)
			self.neighbor_indices.append(md.compute_neighbors(lipids,self.cutoff,[it],periodic=True))
		return None

	def compute_area_per_lipid(self,struc,fr):
		"""
		Compute Solvent Accessible Surface area of membrane WITHIN periodic bounds! 
	
		"""
		from timeit import default_timer as timer
		start = timer()

		sel = "resname DMPC and not type H"
		lipid_ind,lipid_res,lipid_resname = self.get_residues(struc,sel)
		lipids = struc.atom_slice(lipid_ind)
		lipid_ind = [_ for _ in lipid_ind]
		lipid_res = [_ for _ in lipid_res]
		n_lipids = len(set(lipid_res))
		n_atoms_per_lipid = len(lipid_res)/n_lipids

		pbc = np.array(struc.unitcell_lengths)[0]
		surface = []
		area_per_lipid = np.zeros(n_lipids)
	
		# iterate through sel atoms
		for it,atom_i in enumerate(lipid_ind):
			print "atom", it, "of", len(lipid_ind)
			# load neighbors of atom i
			neighbor_indices = self.neighbor_indices[it][0]
			# make mdtraj object containing neighbors
			neighbor_slice = lipids.atom_slice(neighbor_indices)
			# store coordinates and element names of all neighbors
			neighbor_xyz = neighbor_slice.xyz[0] ; neighbor_el = [neighbor_slice.topology.to_dataframe()[0]['element'][i] for i in range(neighbor_slice.n_atoms)]
			# compute radius around atom i (including probe distance)
			r_i = self.r_i[it]
			# define constants
			n_accessible_points = 0
			# iterate through each test point
			for point_k in self.sphere_points:
				is_touching = False 
				point_k = (point_k*r_i + self.coors_i[it][0])[0]
				# cycle through each atom that's near atom_i
				for j in range(len(neighbor_indices)):
					# define probe radius of atom_j		
					r_j = self.probe + ATOMIC_RADII[neighbor_el[j]]
					# measure distance between atom_j and sel1 test point k
					d_jk = np.linalg.norm(self.pbc_dist(neighbor_xyz[j],point_k,pbc))
					# if distance between atom_j and test point is smaller than atom_j radius, then it's NOT accessible
					if d_jk < r_j:
						is_touching = True ; break
				if is_touching == False: 
					n_accessible_points += 1
					surface.append(point_k) 
			area_per_lipid[it//n_atoms_per_lipid] += self.const*n_accessible_points*r_i*r_i
		self.area_per_lipid.append(np.average(area_per_lipid))
		#self.write_pdb(surface,self.outdir+"test_points"+str(fr)+".pdb")
		print "time:", timer()-start, "n points", n_accessible_points
		return None
	
	def MASA_protein(self,struc,sel1="protein",sel2="resname DMPC"):
		"""
		Compute Membrane Accessible Surface Area (MASA) of a protein
		using an adapted Shrake-Rupley algorithm

		Similar to area per lipid calculation, but instead of looking for
		test points that are NOT near other lipid atoms, this looks
		for test points (around protein atoms) that ARE near lipid atoms.
	
		"""
		# constants
		pbc = np.array(struc.unitcell_lengths)[0]
		surface = []
		prot_pts = []
		indlist, reslist, rnlist = self.get_residues(struc,sel1)
	
		# make selections
		sel1_ind = struc.topology.select(sel1) ; sel2_ind = struc.topology.select(sel2)
		sel1_neighbors = md.compute_neighbors(struc,self.cutoff,sel2_ind,sel1_ind)[0]
		sel2_neighbors = md.compute_neighbors(struc,self.cutoff,sel1_ind,sel2_ind)[0]
		areas = np.zeros(len(sel1_ind)) ; first_atom = min(sel1_ind)
	
		# iterate through sel1 atoms that are neighbors of sel2
		for it,atom_i in enumerate(sel1_neighbors):
			# select atom i
			this_atom = struc.atom_slice(struc.topology.select('index '+str(atom_i)))[0]
			# find neighbors of atom i from sel2
			neighbor_indices = md.compute_neighbors(struc,self.cutoff,[atom_i],sel2_ind)[0]
			# make mdtraj object containing neighbors
			neighbor_slice = struc.atom_slice(neighbor_indices)
			# store coordinates and element names of all neighbors
			neighbor_xyz = neighbor_slice.xyz[0] ; neighbor_el = [neighbor_slice.topology.to_dataframe()[0]['element'][i] for i in range(neighbor_slice.n_atoms)]
			# compute radius around atom i (including probe distance)
			r_i = self.probe + ATOMIC_RADII[this_atom.topology.to_dataframe()[0]['element'][0]]
			# define constants
			n_accessible_points = 0 ; is_accessible = True
			# iterate through each sphere point
			for k in range(self.n_sphere_point):
				# move the sphere point and adjust the radius
				point_k = (self.sphere_points[k, :]*r_i + this_atom.xyz[0])[0]
				# cycle through each "sel2" atom that's near sel_1 atom_i
				for j in range(len(neighbor_indices)):
					# define probe radius of sel2 atom_j		
					r_j = self.probe + ATOMIC_RADII[neighbor_el[j]]
					# measure distance between sel2 atom_j and sel1 test point k
					d_jk = np.linalg.norm(self.pbc_dist(neighbor_xyz[j],point_k,pbc))
					# if distance between atom_j and test point is smaller than atom_j radius, then it's NOT accessible
					if d_jk > r_j: 
						is_accessible = False
					else: 
						# save coordinates of protein atom
						prot_pts.append(this_atom.xyz[0][0])
						# save coordinates of surface atom
						surface.append(point_k) 
						# keep track of number of accessible points
						n_accessible_points+=1
						# if the test point is accessible to at least one point, then no need to count it twice
						break
				# until we find an accessible point, keep looking
				if is_accessible == True: break
			#areas[it] = const*n_accessible_points*r_i*r_i
			areas[atom_i-first_atom] = self.const*n_accessible_points*r_i*r_i
		
		ureslist = sorted(set(reslist))
		resareas = np.zeros(len(ureslist))
		first_res = min(ureslist)
	
		for ai,area in enumerate(areas): 
			if area > 0: resareas[reslist[ai]-first_res] += area
		#write_pdb(surface,"test_points.pdb") ; write_pdb(prot_pts,"neigh_points.pdb")
		print "MASA Results!!!!", np.where(np.array(resareas) > 0)[0]
		return np.where(np.array(resareas) > 0)[0]

class lipid_analysis:
	"""
	This class has several things hard-coded, so watch out.

	Calculate the deuterium order parameter and density of lipids using gmx functions.
	You must already have 'gmx order' and 'gmx density' installed for this to work.

	This code doesn't contribute anything new, it's just a simple wrapper
	so that gmx functions can be added into the pipeline.

	"""

	def __init__(self,outdir,name_mod):
		self.outdir = outdir
		self.name_mod = name_mod
		self.first_frame = -1
		self.last_frame = -1
		self.ndx = self.outdir+'index_custom.ndx'
		self.ndx_chain = self.outdir+'chain1.ndx'
		self.order_file = self.outdir+'order'+self.name_mod+'.xvg'
		self.density_file = self.outdir+'density'+self.name_mod+'.xvg'

	def make_ndx(self, tpr):
		if not os.path.isfile(self.ndx):
			print "Creating", self.ndx
			trjconv = subprocess.Popen(['gmx','make_ndx','-f',tpr,'-o',self.ndx],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13 & a C21 | a C22 | a C23 | a C24 | a C25 | a C26 | a C27 | a C28 | a C29 | a C210 | a C211 | a C212 | a C213 | a C214\nname 17 chain1\n13 & a C31 | a C32 | a C33 | a C34 | a C35 | a C36 | a C37 | a C38 | a C39 | a C310 | a C311 | a C312 | a C313 | a C314\nname 18 chain2\nq\n')
			trjconv.wait()
		if not os.path.isfile(self.ndx_chain):
			print "Creating", self.ndx_chain
			trjconv = subprocess.Popen(['gmx','make_ndx','-f',tpr,'-o',self.ndx_chain],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13 & a C21\n 13 & a C22\n 13 & a C23\n 13 & a C24\n 13 & a C25\n 13 & a C26\n 13 & a C27\n 13 & a C28\n 13 & a C29\n 13 & a C210\n 13 & a C211\n 13 & a C212\n 13 & a C213\n 13 & a C214\n\ndel 0-16\nq\n')
			trjconv.wait()
		return None

	def gmx_density(self,xtc, tpr, overwrite,first_frame=-1,last_frame=-1):
		# gmx density -f traj_pbc.xtc -n index.ndx -center
		if not os.path.isfile(self.density_file):
			if first_frame != -1 and last_frame != -1:
				first_frame *= 1000
				last_frame *= 1000
				trjconv = subprocess.Popen(['gmx','density','-f',xtc,'-n',self.ndx, '-s', tpr, '-center', '-o', self.density_file,'-b',str(first_frame),'-e',str(last_frame)],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			else:
				trjconv = subprocess.Popen(['gmx','density','-f',xtc,'-n',self.ndx, '-s', tpr, '-center', '-o', self.density_file],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			trjconv.communicate(b'13\n13\n')
			trjconv.wait()
		return None

	def gmx_order(self,xtc, tpr, overwrite,first_frame=-1,last_frame=-1):
		# gmx order -s topol.tpr -f traj_pbc.xtc -n ac2.ndx -d z -od deuter_ac2.xvg
		if not os.path.isfile(self.order_file):
			if first_frame != -1 and last_frame != -1:
				first_frame *= 1000
				last_frame *= 1000
				trjconv = subprocess.Popen(['gmx','order','-f',xtc,'-n',self.outdir+'chain1.ndx','-nr',self.outdir+'chain1.ndx','-s', tpr,'-d','z','-od', self.order_file,'-b',str(first_frame),'-e',str(last_frame)],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
			else:
				trjconv = subprocess.Popen(['gmx','order','-f',xtc,'-n',self.outdir+'chain1.ndx','-nr',self.outdir+'chain1.ndx','-s', tpr,'-d','z','-od', self.order_file],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
				trjconv.communicate(b'17\n')
				trjconv.wait()
		return None

	def plot_order(self):
		for fil in [self.order_file, self.density_file]:
			x = [] ; y = []
			for line in open(fil):
				try:
					line = line.split()
					line[0]
				except:
					continue
				if line[0] == '#':
					continue
				if len(line) == 2 and '@' not in line[0]:
					x.append(float(line[0]))
					y.append(float(line[1]))

			plt.clf()
			plt.scatter(x,y)
			if fil == self.order_file:
				plt.xlabel('Carbon Number') 
				plt.ylabel('Order Parameter')
				print "plotting", self.outdir+'order_parameter'+self.name_mod+'.png'
				plt.savefig(self.outdir+'order_parameter'+self.name_mod+'.png')
			else:
				plt.xlabel('Z') 
				plt.ylabel('Density')
				print "plotting", self.outdir+'lipid_density'+self.name_mod+'.png'
				plt.savefig(self.outdir+'lipid_density'+self.name_mod+'.png')
		return None

	def lipid_phase_transition(self,struc):
		sel_ind = struc.topology.select('resname DMPC and name 4C31')
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		self.lipid_ids = [i for i in table['resSeq']]
		self.n_lipids = len(self.lipid_ids)
		upper_leaflet = self.lipid_ids[:self.n_lipids/2]
		lower_leaflet = self.lipid_ids[self.n_lipids/2:]
		print upper_leaflet, lower_leaflet
		
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

	def __init__(self,outdir):
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

three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class contacts:
	"""
	Simplest possible method of identifying contacts between proteins and lipids.

	I don't particularly like this method because it's not very accurate yet it's
	still relatively time consuming. The MASA method I use (which may or may not
	be an algorithm other people use and I just don't know it) is much more 
	accurate but even more expensive. 

	Play with the cutoff a bit. I think 0.3 nm is decent, but remember different
	amino acids/different types of interactions have different equilibrium bond
	distances. For example, Arginine makes a closer bonds with lipids through
	charge-charge than tyrosine does through hydrophobic interactions. If your
	cutoff is too small then you'll get lots of arginine contacts and few tyrosine
	contacts.

	"""

	def __init__(self,sel1="protein",sel2="resname DMPC",cutoff=0.3):
		self.sel1 = sel1
		self.sel2 = sel2
		self.cutoff = cutoff
		self.N = -1
		self.membrane_contacts = []
		self.contact_frames = []
		self.first_contact = -1

	def contacts_precalcs(self,struc):
		self.indlist, self.reslist, self.rnlist = self.get_residues(struc)
		a,b,c=self.get_residues(struc,"name CA")
		self.aaseq = [three_to_one[r] for r in c]
		self.sel1_ind = struc.topology.select(self.sel1)
		self.sel2_ind = struc.topology.select(self.sel2)
		self.first_atom = min(self.sel1_ind)
		return None

	def get_residues(self,struc,asel="protein"):
		sel_ind = struc.topology.select(asel)
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		self.N = max(table['resSeq'])
		return table['serial'], table['resSeq'], table['resName']

	def compute_contacts(self,struc):
		sel1_neighbors = md.compute_neighbors(struc,self.cutoff,self.sel2_ind,self.sel1_ind)[0]
		contact_res = sorted(set([self.reslist[i-self.first_atom] for i in sel1_neighbors]))
		for c in contact_res:
			self.membrane_contacts = np.append(self.membrane_contacts,c)
		fr_contact = np.zeros(len(self.aaseq))
		for c in contact_res:
			fr_contact[c-2] = 1
		self.contact_frames.append(fr_contact)
		return None

	def contact_history(self,outdir,name_mod):
		history = np.zeros((np.array(self.contact_frames).shape))
		no_contacts = 0 ; contact_fr = 0
		for fr_,fr in enumerate(self.contact_frames):
			if fr_ == 0: continue
			if any(fr) != 0:
				contact_fr += 1
				# record how many contacts in a row
				inds = np.where(fr != 0)[0]
				inds_prev = np.where(self.contact_frames[fr_-1] != 0)[0]
				for i in inds:
					if i in inds_prev:
						v = history[fr_-1][i]
					else:
						v = 1
					history[fr_][i] += 1 + v
				# first contact
				if self.first_contact == -1:
					self.first_contact = fr_
					print "first contact! Frame:", fr_
					with open(outdir+"contacts_summary" + name_mod + ".txt", 'a') as f:
						f.write("first contact frame: "+str(fr_)+"\n")
			# count frames with no contacts
			else:
				no_contacts += 1
		with open(outdir+"contacts_summary" + name_mod + ".txt", 'a') as f:
			f.write("longest contact lasted: "+str(np.amax(history))+" frames\n")
			f.write("number of frames with contacts: "+str(contact_fr)+"\n")
			f.write("number of frames without contacts :"+str(no_contacts)+"\n\n")
		print "max history", np.amax(history)
		return None

	def hist_restype(self,outdir,name_mod):
		plt.clf()
		restypes = set(self.aaseq)
		hist = np.zeros(len(restypes))
		for mc in self.membrane_contacts:
			resname = self.aaseq[int(mc)-2]
			list_index = list(restypes).index(resname)
			hist[list_index]+=1
		print hist
	
		# zip hist with restypes, sort the *unpacked lists, and then repack them with zip
		hist, restypes = (list(t) for t in zip(*sorted(zip(hist, restypes),reverse=True)))
	
		fig=plt.figure(1)
		ax=fig.add_subplot(111)
		fig.suptitle('Residue/Lipid Contacts')
		ax.set_xlabel('Residue')
		ax.set_ylabel('Number of Occurrences')
	
		ax.set_xticks(np.arange(len(restypes)))
		ax.set_xticklabels(restypes,fontSize=5)
	
		plt.bar(np.arange(len(restypes)),hist,width=1,edgecolor='k')
		plt.xlim(-1,len(restypes))
		plt.savefig(outdir+"contacts_restype" + name_mod + ".png")
		return None

	def plot_contact_hist(self,outdir,name_mod):
		self.contact_history(outdir,name_mod)
		self.hist_restype(outdir,name_mod)
		plt.clf()
		plt.hist(self.membrane_contacts,bins=range(2,self.N+2),edgecolor='black',align='left')
		plt.xticks(range(2,len(self.aaseq)+2), list(self.aaseq))
		print "plotting", outdir+"contacts" + name_mod + ".png"
		plt.savefig(outdir+"contacts" + name_mod + ".png")
		return None
