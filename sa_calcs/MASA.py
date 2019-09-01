from utils import np, md, plt

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

	def __init__(self,outdir,name_mod):
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
	def MASA_get_residues(struc,sel):
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
		lipid_ind,lipid_res,lipid_resname = self.MASA_get_residues(struc,sel)
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
		lipid_ind,lipid_res,lipid_resname = self.MASA_get_residues(struc,sel)
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
		indlist, reslist, rnlist = self.MASA_get_residues(struc,sel1)
	
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

