from lib_handler import np, md, plt, subprocess, os

three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

scale = 2
hg_atoms = '(name N or name C13 or name C14 or name C15 or name C12 or name C11 or name P or name O13 or name O14 or name O11 or name O12 or name C1)'

class mem_analysis:
	def __init__(self,outdir):
		self.area_per_lipid = [] # Area Per Lipid
		self.outdir = outdir
		self.ndx = self.outdir+'index_custom.ndx'
		self.ndx_chain = self.outdir+'chain1.ndx'
		self.av_heights = np.array([])
		self.n_box_counts = np.array([])
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
		else:
			return c

	def heights(self,struc):
		protein_sel = struc.topology.select('protein and backbone')
		prot_coor = np.mean(struc.atom_slice(protein_sel).xyz,axis=1)[0][0:2]
		center_x, center_y = self.snap_to_grid(prot_coor[0]),self.snap_to_grid(prot_coor[1])
		shift_x, shift_y = center_x-self.cx,center_y-self.cy
		for lipid in self.lipid_ids:
			lipid_sel = struc.topology.select(hg_atoms+' and residue ' +str(lipid))
			sel = struc.atom_slice(lipid_sel)
			table, bonds = sel.topology.to_dataframe()
			lc = np.mean(sel.xyz,axis=1)[0]
			x = self.wrap_back(self.snap_to_grid(lc[0])-shift_x,self.nx) 
			y = self.wrap_back(self.snap_to_grid(lc[1])-shift_y,self.ny)
			z = np.mean(sel.xyz.T[2],axis=0)
			self.n_box_counts[x][y] += 1
			self.av_heights[x][y] += z
		return None

	def normalize_heights(self):
		for r,row in enumerate(self.av_heights):
			for c in range(len(row)):
				self.av_heights[r][c]/=float(self.n_box_counts[r][c])
		print "av_heights!", self.av_heights
		fig, ax = plt.subplots()
		cax = ax.imshow(self.av_heights)
		fig.colorbar(cax)
		plt.savefig(self.outdir+'average_heights_centered.png')
		return None

	def heights_precalcs(self,struc):
		sel_ind = struc.topology.select('resname DMPC and name C21')
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		self.lipid_ids = [i for i in table['resSeq']]
		self.n_lipids = len(self.lipid_ids)
		self.lipid_ids = self.lipid_ids[:self.n_lipids/2] # upper leaflet only
		self.nx = scale*int(np.ceil(struc.unitcell_lengths[0][0]))
		self.ny = scale*int(np.ceil(struc.unitcell_lengths[0][1]))
		self.av_heights = np.zeros((self.nx,self.ny))
		self.n_box_counts = np.zeros((self.nx,self.ny))
		self.cx = self.nx/2
		self.cy = self.ny/2
		return None

	def compute_area_per_lipid(self,struc):
		self.area_per_lipid.append(struc.unitcell_lengths[0][0]*struc.unitcell_lengths[0][1]/77.0)
		return None

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

	def gmx_density(self,xtc, tpr):
		# gmx density -f traj_pbc.xtc -n index.ndx -center
		print "Computing density for all DMPC"
		trjconv = subprocess.Popen(['gmx','density','-f',xtc,'-n',self.ndx, '-s', tpr, '-center', '-o', self.outdir+'density.xvg'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		trjconv.communicate(b'13\n13\n')
		trjconv.wait()

	def gmx_order(self,xtc, tpr):
		# gmx order -s topol.tpr -f traj_pbc.xtc -n ac2.ndx -d z -od deuter_ac2.xvg
		print "Computing order for all chain1"
		trjconv = subprocess.Popen(['gmx','order','-f',xtc,'-n',self.outdir+'chain1.ndx','-nr',self.outdir+'chain1.ndx','-s', tpr,'-d','z','-od', self.outdir+'order.xvg'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		trjconv.communicate(b'17\n')
		trjconv.wait()

class contacts:

	def __init__(self,sel1="protein",sel2="resname DMPC",cutoff=0.23):
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

	def plot_contact_hist(self,outdir,name_mod):
		self.contact_history(outdir,name_mod)
		self.hist_restype(outdir,name_mod)
		plt.clf()
		plt.hist(self.membrane_contacts,bins=range(2,self.N+2),edgecolor='black',align='left')
		plt.xticks(range(2,len(self.aaseq)+2), list(self.aaseq))
		print "plotting", outdir+"contacts" + name_mod + ".png"
		plt.savefig(outdir+"contacts" + name_mod + ".png")

