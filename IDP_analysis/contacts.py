from lib_handler import np, md, plt

"""
Merged with sa_mem

"""

class contacts:

	def __init__(self,sel1="protein",sel2="resname DMPC",cutoff=0.3):
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
		print "!!!", self.membrane_contacts
		plt.hist(self.membrane_contacts,edgecolor='black',align='left')
		#plt.hist(self.membrane_contacts,bins=range(1,self.N+1),edgecolor='black',align='left')
		plt.savefig(outdir+"contacts" + name_mod + ".png")

