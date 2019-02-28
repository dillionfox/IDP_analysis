from lib_handler import np, md, plt

class mem_analysis:
	def __init__(self):
		self.APL = [] # Area Per Lipid

	def area_per_lipid(self,struc):
		self.APL.append(struc.unitcell_lengths[0][0]*struc.unitcell_lengths[0][1]/77.0)
		return None

class contacts:

	def __init__(self,sel1="protein and backbone",sel2="resname DMPC",cutoff=0.2):
		self.sel1 = sel1
		self.sel2 = sel2
		self.cutoff = cutoff
		self.N = -1
		self.membrane_contacts = []

	def compute_contacts(self,struc):
		indlist, reslist, rnlist = self.get_residues(struc,self.sel1)
		sel1_ind = struc.topology.select(self.sel1)
		sel2_ind = struc.topology.select(self.sel2)
		sel1_neighbors = md.compute_neighbors(struc,self.cutoff,sel2_ind,sel1_ind)[0]
		first_atom = min(sel1_ind)
		contact_res = sorted(set([reslist[i-first_atom] for i in sel1_neighbors]))
		for c in contact_res:
			self.membrane_contacts = np.append(self.membrane_contacts,c)
		return None

	def get_residues(self,struc,sel):
		sel_ind = struc.topology.select(sel)
		sel = struc.atom_slice(sel_ind)
		table, bonds = sel.topology.to_dataframe()
		self.N = max(table['resSeq'])
		return table['serial'], table['resSeq'], table['resName']

	def plot_contact_hist(self,outdir,name_mod):
		plt.clf()
		plt.hist(self.membrane_contacts,edgecolor='black',align='left')
		plt.hist(self.membrane_contacts,bins=range(0,self.N+2),edgecolor='black',align='left')
		print "plotting", outdir+"contacts" + name_mod + ".png"
		plt.savefig(outdir+"contacts" + name_mod + ".png")

