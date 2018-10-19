#!/usr/bin/python
import numpy as np
#import sys
#import MDAnalysis
import mdtraj as md
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt

### Protein/Lipid SASA ###
global ATOMIC_RADII
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

def plot_masa(l,sequence,name):
	plt.clf()
	flat_list = [item for sublist in l for item in sublist]

	###

	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	fig.suptitle('Residue/Lipid Contacts')
	ax.set_xlabel('Residue Number')
	ax.set_ylabel('Number of Occurrences')

	ax.set_xticks(np.arange(len(sequence)))
	ax.set_xticklabels(sequence,fontSize=5)

	plt.xlim(-1,len(sequence))

	###

	#plt.hist(flat_list,align='mid',rwidth=5)
	hist, bins = np.histogram(flat_list, bins = range(0,len(sequence)))
	plt.bar(bins[:-1], hist, width=1, edgecolor='k')
	
	plt.savefig(name+"_contacts.png")
	return None

def generate_sphere_points(n):
	"""
	generate coordinates on a sphere using the Golden Section Spiral algorithm

	"""
	points = np.zeros((n, 3))
	inc = np.pi * (3 - np.sqrt(5))
	offset = 2/float(n)
	for k in range(int(n)):
		y = k*offset - 1 + (offset/2)
		r = np.sqrt(1 - y*y)
		phi = k*inc
		points[k, :] = [np.cos(phi) * r, y, np.sin(phi) * r]
	return points

def pbc_dist(r_i,r_j,pbc):
	d = r_i - r_j
	for q in range(3):
		if abs(d[q]) > pbc[q]/2.0:
			if d[q] > 0: d[q] -= pbc[q]
			else: d[q] += pbc[q]
	return d

def get_residues(struc,sel):
	sel_ind = struc.topology.select(sel)
	sel = struc.atom_slice(sel_ind)
	table, bonds = sel.topology.to_dataframe()
	return table['serial'], table['resSeq'], table['resName']

def MASA(struc,sel1 = "protein", sel2 = "resname DOPC or resname DOPS or resname POPC"):
	"""
	Compute Membrane Accessible Surface Area (MASA) of a protein
	using an adapted Shrake-Rupley algorithm
	sel1 = "protein and not (name =~ 'H.*')"

	"""
	# constants
	n_sphere_point = 960 ; cutoff = 0.4 ; probe = 0 ; sphere_points = generate_sphere_points(n_sphere_point)
	const = 4.0*np.pi/n_sphere_point ; surface = [] ; pbc = np.array(struc.unitcell_lengths)[0] ; prot_pts = []
	indlist, reslist, rnlist = get_residues(struc,sel1)

	# make selections
	sel1_ind = struc.topology.select(sel1) ; sel2_ind = struc.topology.select(sel2)
	sel1_ind = struc.topology.select(sel1) ; sel2_ind = struc.topology.select(sel2)
	sel1_neighbors = md.compute_neighbors(struc,cutoff,sel2_ind,sel1_ind)[0]
	sel2_neighbors = md.compute_neighbors(struc,cutoff,sel1_ind,sel2_ind)[0]
	areas = np.zeros(len(sel1_ind)) ; first_atom = min(sel1_ind)

	# iterate through sel1 atoms that are neighbors of sel2
	for it,atom_i in enumerate(sel1_neighbors):
		# select atom i
		this_atom = struc.atom_slice(struc.topology.select('index '+str(atom_i)))[0]
		# find neighbors of atom i from sel2
		neighbor_indices = md.compute_neighbors(struc,cutoff,[atom_i],sel2_ind)[0]
		# make mdtraj object containing neighbors
		neighbor_slice = struc.atom_slice(neighbor_indices)
		# store coordinates and element names of all neighbors
		neighbor_xyz = neighbor_slice.xyz[0] ; neighbor_el = [neighbor_slice.topology.to_dataframe()[0]['element'][i] for i in range(neighbor_slice.n_atoms)]
		# compute radius around atom i (including probe distance)
		r_i = probe + ATOMIC_RADII[this_atom.topology.to_dataframe()[0]['element'][0]]
		# define constants
		n_accessible_points = 0 ; is_accessible = True
		# iterate through each sphere point
		for k in range(n_sphere_point):
			# move the sphere point and adjust the radius
			point_k = (sphere_points[k, :]*r_i + this_atom.xyz[0])[0]
			# cycle through each "sel2" atom that's near sel_1 atom_i
			for j in range(len(neighbor_indices)):
				# define probe radius of sel2 atom_j		
				r_j = probe + ATOMIC_RADII[neighbor_el[j]]
				# measure distance between sel2 atom_j and sel1 test point k
				d_jk = np.linalg.norm(pbc_dist(neighbor_xyz[j],point_k,pbc))
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
		areas[atom_i-first_atom] = const*n_accessible_points*r_i*r_i
	ureslist = sorted(set(reslist)) ; resareas = np.zeros(len(ureslist)) ; first_res = min(ureslist)
	for ai,area in enumerate(areas): 
		if area > 0: resareas[reslist[ai]-first_res] += area
	#write_pdb(surface,"test_points.pdb") ; write_pdb(prot_pts,"neigh_points.pdb")
	return np.where(np.array(resareas) > 0)[0]

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


### interface probe ###
vec_to_int = {'x':0, 'y':1, 'z':2}
int_to_vec = {0:'x', 1:'y', 2:'z'}
interface_vector = 'z'

def probe(fr,GRO,XTC,cutoff,probe_radius):
	print "frame:", fr, "cutoff:", cutoff
	uni = MDAnalysis.Universe(GRO,XTC)
	uni.trajectory[fr]
	bilayer = uni.select_atoms('resname POPC or resname DOPC or resname DOPS')
	selstring = 'around ' + str(cutoff) +  ' global (name BB)'
	lipids = bilayer.select_atoms(selstring)
	lipid_av = lipids.positions[:,vec_to_int[interface_vector]].mean(axis=0)

	protein_all = uni.select_atoms('protein ')
	resids = np.zeros(len(protein_all.residues.resids))
	selstring = 'around ' + str(cutoff) +  ' global (resname POPC or resname DOPC or resname DOPS)'
	protein = protein_all.select_atoms(selstring)
	protein_av = protein.positions[:,vec_to_int[interface_vector]].mean(axis=0)

	if protein_av > lipid_av: 	# if protein is ABOVE bilayer
		ineq = '<'		# then we are looking for water BELOW protein
		ineq_ = '>'
	else: 				# if protein is BELOW bilayer
		ineq = '>'		# then we are looking for water ABOVE protein
		ineq_ = '<'

	for i in protein.residues.resids:
		pc = protein_all.positions[i-1]
		iv = interface_vector # shorter name
		dl = probe_radius # shorter name
		v0 = vec_to_int[iv]; v1 = (v0+1)%3; v2 = (v0+2)%3
		iv1 = int_to_vec[v1]; iv2 = int_to_vec[v2]

		selstring = '(resname W or name BB) and (prop '+str(iv)+ineq+str(pc[v0])
		selstring += ' and prop ' + str(iv) + ineq_ + str(lipid_av) 
		selstring += ' and prop ' + str(iv1) + ' > ' + str(pc[v1]-dl) 
		selstring += ' and prop ' + str(iv1) + ' < ' + str(pc[v1]+dl) 
		selstring += ' and prop ' + str(iv2) + ' > ' + str(pc[v2]-dl) 
		selstring += ' and prop ' + str(iv2) + ' < ' + str(pc[v2]+dl) 
		selstring += ')'
		water_sel = uni.select_atoms(selstring)

		if len(water_sel) == 0:
			print "CONTACT:", i
			resids[i-1] = 1
	return resids 

def probe_make_map(dat, sequence):
	from matplotlib.ticker import MultipleLocator

	def average_set(curr_list):
		cutoff = 5
	        hist=[]
	        unique_list=sorted(set(curr_list))
	        for el in unique_list:
	                occ=curr_list.count(el)
	                if occ > cutoff:
	                        hist.append(el)
		return hist
	
	fontSize = 10
	res_per_frame = []	
	seq = []
	curr_list = []
	frame_bin = []
	fr_count = 0

	tdat = np.transpose(dat)

	for frame in tdat:
		residues = np.where(frame!=0)[0]
	        for r in residues:
			curr_list.append(r)
	        if fr_count%25==0:
			print "HERE", curr_list
	                hist=average_set(curr_list)
	                for i in range(len(hist)):
	                        frame_bin.append(fr_count)
	                        seq.append(hist[i])
	                res_per_frame.append(len(hist))
	                curr_list=[]
	        fr_count+=1

        fig=plt.figure(1)
        plot=fig.add_subplot(111)
        plt.scatter(frame_bin,seq,s=30,marker='o',lw=0.0)
        plot.tick_params(axis='both', which='minor', labelsize=fontSize)
        #plot.yaxis.set_minor_locator(MultipleLocator(1))
        #plt.grid(b=True,which='major',linestyle='-')
        #plt.grid(b=True,which='minor',linestyle='--',alpha=0.3)
        plot.set_xlabel('Time (ns)',size=fontSize)
        plot.set_ylabel('Residue',size=fontSize)
        plot.set_yticks(np.arange(0,len(sequence)))
        plot.set_yticklabels(sequence,fontSize=fontSize)
        #plt.ylim(26.5,35.5)   
        #plt.xlim(0,900)


def probe_make_hist_by_restype(x, sequence):
	restypes = set(sequence)

	hist = np.zeros(len(restypes))
	resnum = 0
	for i in x:
		resname = sequence[resnum]
		list_index = list(restypes).index(resname)
		hist[list_index]+=sum(i)/sequence.count(resname)
		print resname, sequence.count(resname)
		resnum+=1
	print hist

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

def probe_make_hist(x, sequence, XTC):
	plt.clf()
	m,n=x.shape
	hist = np.zeros(m)
	resnum = 0
	for i in x:
		hist[resnum]=sum(i)
		resnum+=1

	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	fig.suptitle('Residue/Lipid Contacts')
	ax.set_xlabel('Residue Number')
	ax.set_ylabel('Number of Occurrences')

	ax.set_xticks(np.arange(len(sequence)))
	ax.set_xticklabels(sequence,fontSize=5)

	plt.bar(np.arange(m),hist,width=1,edgecolor='k')
	plt.xlim(-1,len(sequence))
	plt.savefig(XTC.split('.')[0]+"_res_hist.png")

def interface_probe(GRO,XTC,skip_frames,first_frame,last_frame,nthreads,cutoff,probe_radius,sequence):
	print "loading...", XTC
	uni = MDAnalysis.Universe(GRO,XTC)
	print "complete" 
	if last_frame == 'last':
		last_frame = len(uni.trajectory)
	protein = uni.select_atoms('protein')
	nres = len(protein.residues)
	first_resid = protein.residues.resids[0]
	uni = []

	frames = range(first_frame,last_frame,skip_frames)
	resids_per_frame = np.zeros((len(frames),nres))
	contact_frames = np.zeros(len(frames))
	if nthreads == 1:
		resids_per_frame = [probe(fr,GRO,XTC,cutoff,probe_radius) for fr in frames]
	elif nthreads > 1:
		from joblib import Parallel,delayed
		from joblib.pool import has_shareable_memory
		#resids_per_frame = Parallel(n_jobs=nthreads)(delayed(probe,has_shareable_memory)(fr) for fr in frames)
	else:
		print "nthreads can't equal", nthreads

	#print resids_per_frame
	resids_per_frame = np.transpose(resids_per_frame)
	probe_make_hist(resids_per_frame, sequence, XTC)
	#probe_make_map(resids_per_frame, sequence)
	#probe_make_hist_by_restype(resids_per_frame, sequence)

if __name__ == "__main__":

	# interface probe
	GRO = sys.argv[1]
	XTC = sys.argv[2]
	skip_frames = 1
	first_frame = 0
	last_frame = 'last'
	nthreads = 12
	cutoff = 5
	probe_radius = 1.0
	sequence = ['R','Q','M','N','F','P','E','R','S','M','D','M','S','N','L','Q','M','D','M','Q','G','R','W','M','D','M','Q','G','R','Y','T','N','P','F','N']
	interface_probe(GRO,XTC,skip_frames,first_frame,last_frame,nthreads,cutoff,probe_radius,sequence)
