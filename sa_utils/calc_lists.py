"""
This object keeps track of all known calculations and their dependencies

"""

class calc_lists:

	def __init__(self):
	#---Every known calculation
		self.known_calcs = ['Rg', 'SASA', 'EED', 'Asph', 'rama', 'cmaps', 'PCA', 'gcmaps',\
				'XRg','SS', 'chain', 'score','flory', 'centroids', 'Gyr', \
				'surface_contacts', 'rmsd', 'probe', 'MASA', 'calibur', \
				'diffusion', 'contact_types', 'contact_residues', 'membrane_contacts',\
				'area_per_lipid','persistence_length','membrane_analysis','av_heights',\
				'interdigitation']

		#---Calculations that depend on each other
		self.deps = {'Rg':['Gyr'], 'Asph':['Gyr'], 'surface_contacts':['SASA'],\
				'contact_residues':['cmaps'], 'contact_types':['cmaps'],\
				'PCA':['Gyr','Rg','SASA','EED','Asph'],'chain':['Gyr','Rg','EED'],\
				}

		#---Every known post-analysis function
		self.post_analysis = ['Rg', 'cmaps', 'gcmaps', 'surface_contacts', 'SS' , 'PCA',\
				'flory', 'chain', 'centroids', 'rama', 'MASA', 'diffusion',\
				 'contact_residues', 'contact_types', 'membrane_contacts',\
				'membrane_analysis','av_heights']

		self.traj_master_list = ['rmsd','calibur','probe','persistence_length','membrane_analysis',\
				'interdigitation']

		self.precalcs_master_list = ['membrane_contacts','membrane_analysis','av_heights','area_per_lipid',\
				'interdigitation']
