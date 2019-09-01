from utils import np, md, plt

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

	def __init__(self,outdir,name_mod,sel1="protein",sel2="resname DMPC",cutoff=0.32):
		self.sel1 = sel1
		self.sel2 = sel2
		self.cutoff = cutoff
		self.N = -1
		self.membrane_contacts = []
		self.contact_frames = []
		self.first_contact = -1

	def contacts_precalcs(self,struc):
		self.indlist, self.reslist, self.rnlist = self.get_residues(struc,"protein")
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
		try:
			self.contact_history(outdir,name_mod)
		except:
			print "can't load contact history"
		self.hist_restype(outdir,name_mod)
		plt.clf()
		plt.hist(self.membrane_contacts,bins=range(2,self.N+2),edgecolor='black',align='left')
		plt.xticks(range(2,len(self.aaseq)+2), list(self.aaseq))
		print "plotting", outdir+"contacts" + name_mod + ".png"
		plt.savefig(outdir+"contacts" + name_mod + ".png")
		return None

