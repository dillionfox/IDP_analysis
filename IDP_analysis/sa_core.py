import numpy as np
import mdtraj as md
import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15}
mpl.rc('font', **font)
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

"""
gyration_tensor(coors): L
compute_Rg(L): Rg
compute_Asph(L): Asph
av_SS(SS): None
plot_SS(H,E,C,nframes,nres,outdir,name_mod): None
contact_maps(coors): dist_array
gremlin_contact_maps(dist): dist_array
surface_contacts(struc,SASA): cmap
av_cmaps(cmaps,nres,resnames,outdir,name_mod,mtype="NULL"): None
compute_XRg(PDB): Rg
scatterplot(X,Y,xlabel,ylabel,filename,outdir,name_mod): None

"""

global hydrophobic ; hydrophobic = np.array(['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'])
global hydrophilic ; hydrophilic = np.array(['SER', 'THR', 'ASN', 'GLN'])
global poscharge   ; poscharge =   np.array(['ARG', 'LYS'])
global negcharge   ; negcharge =   np.array(['ASP', 'GLU', 'HIS'])
global aromatic    ; aromatic =    np.array(['TYR','PHE','TRP'])
global methionine  ; methionine =  np.array(['MET'])
global aliphatic   ; aliphatic =   np.array(['GLY', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'CYS'])
global aa          ; aa =          np.array(['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'ASN', 'GLN', 'ARG', 'LYS', 'ASP', 'GLU', 'HIS', 'PRO', 'GLY', 'CYS'])

class sa_core:
	def __init__(self):
		self.EED = [] 			# end-to-end distance
		self.Rg = [] 			# radius of gyration
		self.XRg = [] 			# predicted x-ray radius of gyration (from crysol)
		self.SASA = [] 			# solvent accessible surface area
		self.Asph = [] 			# asphericity 
		self.SS = []			# secondary structure (per residue)
		self.cmaps = []			# contact maps
		self.gcmaps = []		# contact maps for specific residues
		self.scmaps = []		# contact maps for residues on surface only
		self.seq = None

	def get_seq(self,struc):
		"""
		Extract sequence from input. MDTraj makes this a little annoying

		"""
		reslist = struc.atom_slice(struc.topology.select("protein and name CA"))[0].topology.to_dataframe()[0]['resName']
		self.seq = [reslist[i] for i in range(len(reslist))]
	
	@staticmethod
	def gyration_tensor(coors):
		"""
		3x3 tensor. It's just an outerproduct of r and r. COORDINATES MUST BE CENETERED/
		S_nm = 1/N Sum_i { r(i)_n r(i)_m }
	
		"""
	        coors -= np.mean(coors,axis=0)
	        S = np.sum([np.einsum('i,j->ij',c,c) for c in coors], axis=0)/(len(coors))
	        L = np.linalg.eig(S)[0]
		return np.sqrt(L)
	
	def compute_Rg(self,L):
		"""
		L's are sqrt(eigenvalues of gyration tensor)
		Rg^2 = L1^2 + L2^2 + L3^2
	
		"""
		self.Rg.append(np.sqrt(np.dot(L,L)))
	
	def compute_Asph(self,L):
	        """
		Compute the Asphericitiy
		Eq From: Simulation Analysis of the Temperature Dependence of Lignin Structure and Dynamics
	
	        """
		self.Asph.append(((L[0]-L[1])**2+(L[1]-L[2])**2+(L[0]-L[2])**2)/(2*sum(L)**2))
		#return ((L[0]-L[1])**2+(L[1]-L[2])**2+(L[0]-L[2])**2)/(2*sum(L)**2)
	
	def av_SS(self,SS,outdir,name_mod):
		"""
		Average over all Seconday Structure
	
		"""
	
		nframes = len(SS) ; nres = len(SS[0][0])
		H = np.zeros((nres,nframes)) # alpha helix
		E = np.zeros((nres,nframes)) # extended beta sheet
		C = np.zeros((nres,nframes)) # unstructured coil
	
		fr = 0
		for ss in SS:
			h = np.where(ss[0] == 'H')[0]
			for hh in h:
				H[hh][fr] = 1
	
			e = np.where(ss[0] == 'E')[0]
			for ee in e:
				E[ee][fr] = 1
	
			c = np.where(ss[0] == 'C')[0]
			for cc in c:
				C[cc][fr] = 1
			fr+=1
	
		np.savetxt(outdir+"SS_H" + name_mod + "_raw.npy",H)
		np.savetxt(outdir+"SS_E" + name_mod + "_raw.npy",E)
		np.savetxt(outdir+"SS_C" + name_mod + "_raw.npy",C)
	
		#super(sa_core, self).plot_SS(H,E,C,nframes,nres,outdir,name_mod)
		self.plot_SS(H,E,C,nframes,nres,outdir,name_mod)
		return None

	@staticmethod
	def plot_SS(H,E,C,nframes,nres,outdir,name_mod):
		"""
		Plot average secondary structure per residue with error bars
	
		"""
		H_av = np.zeros(nres) ; H_err = np.zeros(nres) 
		E_av = np.zeros(nres) ; E_err = np.zeros(nres) 
		C_av = np.zeros(nres) ; C_err = np.zeros(nres) 
		for ind in range(nres):
			H_av[ind] = sum(H[ind,:])/nframes ; H_err[ind] = np.std(H[ind,:])
			E_av[ind] = sum(E[ind,:])/nframes ; E_err[ind] = np.std(E[ind,:])
			C_av[ind] = sum(C[ind,:])/nframes ; C_err[ind] = np.std(C[ind,:])
	
		plt.clf()
		plt.plot(H_av,label='Helix',c='r') ;      plt.errorbar(range(nres), H_av,yerr=H_err/(np.sqrt(nframes)-1), fmt='o',color='r')
		plt.plot(E_av,label='Beta Sheet',c='b') ; plt.errorbar(range(nres), E_av,yerr=E_err/(np.sqrt(nframes)-1), fmt='o',color='b')
		plt.plot(C_av,label='Coil',c='k') ;       plt.errorbar(range(nres), C_av,yerr=C_err/(np.sqrt(nframes)-1), fmt='o',color='k')
		plt.legend(loc=1)
		plt.savefig(outdir+"SS" + name_mod + ".png")
		np.savetxt(outdir+"SS_H_av" + name_mod + ".npy",[H_av,H_err])
		np.savetxt(outdir+"SS_E_av" + name_mod + ".npy",[E_av,E_err])
		np.savetxt(outdir+"SS_C_av" + name_mod + ".npy",[C_av,C_err])
		return None
	
	@staticmethod
	def contact_maps(coors):
		"""
		MDAnalysis takes all of the fun out of this one, but the mdad routines are all done within a KD Tree so
		it's faster than anything I would write.
	
		"""
		import MDAnalysis.analysis.distances as mdad
		return mdad.distance_array(coors, coors)
	
	@staticmethod
	def contact_residues(dist,seq,nframes):
		"""
		Tertiary contacts by residue type
	
		try other methods: 	positive-negative charge interactions
					aromatic-methionine
	
		"""
		seq = np.array(seq) ; dist = np.array(dist) ; inds = []
		ut = np.triu_indices(dist.shape[0]) ; dist[ut] = dist.T[ut] # symmetrize matrix. diagonal already 0.
		# iterate through residue
		for r1 in aa:
			for r2 in aa:
				# extract indices of each atom for each type. Complicated because 'flatten' doesnt seem to work
				ind1 = np.where(seq == r1)[0]
				ind2 = np.where(seq == r2)[0]
				if np.sum(dist[ind1].T[ind2]) > 0.9:
					print r1, r2, np.sum(dist[ind1].T[ind2]), len(ind1)*len(ind2), np.sum(dist[ind1].T[ind2])/(len(ind1)*len(ind2))
				#print ind1, ind2, "\n", "\n"
		exit()
		return None
	
	@staticmethod
	def contact_types(dist,seq,nframes):
		"""
		Tertiary contacts by residue type
	
		try other methods: 	positive-negative charge interactions
					aromatic-methionine
	
		"""
		seq = np.array(seq) ; dist = np.array(dist) ; inds = []
		ut = np.triu_indices(dist.shape[0]) ; dist[ut] = dist.T[ut] # symmetrize matrix. diagonal already 0.
		rtype1 = [poscharge, methionine,  aromatic,  aromatic]
		rtype2 = [negcharge,   aromatic,  aromatic, aliphatic]
		labels = [    '+ -',  'met-aro', 'aro-aro', 'aro-ali']
		#rtype1 = [negcharge, poscharge, methionine,   aromatic]
		#rtype2 = [poscharge, negcharge,   aromatic, methionine]
		#labels = [    '- +',     '+ -',  'met-aro',  'aro-met']
		# iterate through residue type
		for r1,r2,l in zip(rtype1,rtype2,labels):
			# extract indices of each atom for each type. Complicated because 'flatten' doesnt seem to work
			ind1 = [item for sublist in [np.where(seq == a)[0] for a in r1] for item in sublist]
			ind2 = [item for sublist in [np.where(seq == a)[0] for a in r2] for item in sublist]
			print l, np.sum(dist[ind1].T[ind2]), len(ind1)*len(ind2), np.sum(dist[ind1].T[ind2])/(len(ind1)*len(ind2))
			#print ind1, ind2, "\n", "\n"
		return None
	
	def gremlin_contact_maps(self,dist):
		"""
		Specifically search for contacts that were used as constraints derived by GREMLIN
	
		"""
	
		print dist
		contact_cutoff = 10
		gremlin = [[6,13],[9,22],[15,19],[14,18],[3,11],[34,40],[3,23],[36,40],[9,13],[25,28],[12,15],[11,23], \
			[26,35],[12,18],[2,5],[17,21],[14,22],[6,9],[41,44],[15,18],[25,30],[9,16],[29,32],[30,33],[6,16]] 
	
		contacts = np.zeros(dist.shape)
		for n in range(dist.shape[0]):
			for m in range(dist.shape[1]):
				if dist[n][m] < contact_cutoff and (gremlin.count([n,m]) == 1 or gremlin.count([m,n]) == 1):
					contacts[n][m] = 1
		self.gcmaps.append(contacts)
	
	@staticmethod
	def surface_contacts(struc,SASA):
		"""
		SASA per residue from SASA per atom
	
		"""
		plot_single = False
	
		cmap_cutoff = 5.0 # nm
		sasa_cutoff = 0.5
		# identify surface residues
		S = []
		N = np.array(range(struc.n_residues))
		for r in N:
			ind = struc.topology.select('resid '+str(r))
			SASA_per_res = np.sum(SASA[0][ind[0]:ind[-1]],axis=0)
			if SASA_per_res > sasa_cutoff:
				S.append(SASA_per_res)
			else:
				S.append(0.0)
		#import ipdb; ipdb.set_trace()
		surf_res = N[np.where(np.array(S) > 0)]
		cmap = contact_maps(struc.atom_slice(struc.topology.select('name CA'))[0].xyz[0])
		for n in N:
			for m in N:
				if (n not in surf_res or m not in surf_res) or n < m:
					cmap[n][m] = 0
				elif cmap[n][m] < cmap_cutoff:
					cmap[n][m] = 1
		if plot_single == True:
			for it,rn in enumerate(resnames):
				if rn in hydrophobic:
					plt.axhline(y=it,c='grey',linewidth=2)
					plt.axvline(x=it,c='grey',linewidth=2)
				elif rn in hydrophilic:
					plt.axhline(y=it,c='g',linewidth=     2)
					plt.axvline(x=it,c='g',linewidth=     2)
				elif rn in poscharge:
					plt.axhline(y=it,c='b',linewidth=     2)
					plt.axvline(x=it,c='b',linewidth=     2)
				elif rn in negcharge:
					plt.axhline(y=it,c='r',linewidth=     2)
					plt.axvline(x=it,c='r',linewidth=     2)
				else:
					print "unknown restype:", rn
			scat = np.array(zip(np.where(np.array(cmap) > 0)[0], np.where(np.array(cmap) > 0)[1])).T
			plt.scatter(scat[0], scat[1],marker='s',c='k',s=90)
			#plt.imshow(cmap,cmap='Greys')
			plt.show()
		return cmap
	
	@staticmethod
	def av_cmaps(cmaps,nres,resnames,outdir,name_mod,mtype="NULL"):
		"""
		Average contact maps. This function has the contact cutoff hard-coded. This
		function works for GREMLIN contact maps too.
	
		"""
		plt.clf()
		nframes = len(cmaps) 
		av = np.zeros((cmaps[0].shape))
	
		# save cmaps to npy file. Data must first be reshaped.
		if mtype == "NULL":
			cmaps = np.array(cmaps)
			resh = cmaps.reshape(cmaps.shape[0],cmaps.shape[1]*cmaps.shape[2])
			np.savetxt(outdir+"CMAPS" + name_mod + "_raw.npy",resh)
	
			for i in range(cmaps[0].shape[0]):
				for j in range(cmaps[0].shape[1]):
					# for each element of the matrix
					if j > i: # don't compute things twice
						l = []
						for c in cmaps:
							# for each map, determine if there was a contact at that position
							if c[i][j] < 0.7: # nm
								l.append(1)
							else:
								l.append(0)
						av[i][j] = np.std(l)/(np.sqrt(nframes)-1)
						av[j][i] = np.mean(l)
					# dont consider contacts from neighbors
					if i == j or abs(i-j) <= 2:
						av[i][j] = 0
						av[j][i] = 0
		else:
			for m in range(nres):
				for n in range(nres):
					for fr in range(nframes):
						av[n][m] += cmaps[fr][m][n]
			av/=nframes
		fig, ax = plt.subplots()
		plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
		if mtype == "gremlin":
			im = ax.imshow(av, cmap='PuBu')
			cbar = fig.colorbar(im)
			ax.set_title("Average Contact Maps from Rosetta+Gremlin Output")
			plt.savefig(outdir+"gremlin_compare_CMAPS" + name_mod + ".png")
			np.savetxt(outdir+"gremlin_CMAPS" + name_mod + "_av.npy",av)
		elif mtype == "surface":
			hydrophobic = ['GLY', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'CYS']
			hydrophilic = ['SER', 'THR', 'ASN', 'GLN', 'HIS']
			poscharge =   ['ARG', 'LYS']
			negcharge =   ['ASP', 'GLU']
	
			for it,rn in enumerate(resnames):
				if rn in hydrophobic:
					plt.axhline(y=it,c='yellow',linewidth=1.5)
					plt.axvline(x=it,c='yellow',linewidth=1.5)
				elif rn in hydrophilic:
					plt.axhline(y=it,c='g',linewidth=1.5)
					plt.axvline(x=it,c='g',linewidth=1.5)
				elif rn in poscharge:
					plt.axhline(y=it,c='b',linewidth=1.5)
					plt.axvline(x=it,c='b',linewidth=1.5)
				elif rn in negcharge:
					plt.axhline(y=it,c='r',linewidth=1.5)
					plt.axvline(x=it,c='r',linewidth=1.5)
				else:
					print "unknown restype:", rn
			ax.set_title("Average Contact Maps of Surface Residues")
			im = ax.imshow(av, cmap='Greys')
			cbar = fig.colorbar(im)
			plt.savefig(outdir+"surface_CMAPS" + name_mod + ".png")
			np.savetxt(outdir+"surface_CMAPS" + name_mod + "_av.npy",av)
		else:
			im = ax.imshow(av)
			cbar = fig.colorbar(im)
			ax.set_title("Average Contact Maps")
			plt.savefig(outdir+"CMAPS" + name_mod + ".png")
			np.savetxt(outdir+"CMAPS" + name_mod + "_av.npy",av)
		return av
	
	@staticmethod
	def compute_XRg(PDB):
		"""
		X-Ray experiments return higher values of Rg because they include some of the water in the shell. The EMBL
		Program "Crysol" computes a theoretical scattering curve for a protein and returns the Rg.
	
		"""
		import subprocess
	
		f = PDB.split('.')[0]
		FNULL = open(os.devnull, 'w')
		subprocess.call(['crysol',f+'.pdb'], stdout=FNULL, stderr=subprocess.STDOUT)
		for line in open(f+'00.log'):
			if "Rg ( Atoms - Excluded volume + Shell ) ................. :" in line:
				os.remove(f+'00.log') ; os.remove(f+'00.alm') ; os.remove(f+'00.int')
				return float(line.split(' : ')[1])
	
	@staticmethod
	def scatterplot(X,Y,xlabel,ylabel,filename,outdir,name_mod):
		"""
		General scatter plot function
	
		"""
		plt.clf()
		plt.scatter(X,Y)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.savefig(outdir+filename + name_mod + ".png")
		np.savetxt(outdir+ xlabel + name_mod + ".npy",x)
		np.savetxt(outdir+ ylabel + name_mod + ".npy",y)
		return None
